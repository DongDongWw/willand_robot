# Trajectory Tracking using MPC of Ackermann Model

## 1 建模
### 1.1 控制变量 & 状态变量
$$
\begin{align*}
& \mathbf{x} = (x, y, \theta, v) \\
& \mathbf{u} = (\omega, a)
\end{align*}
$$
### 1.2 MPC范式
- 代价函数：一般代价和终端代价，一般代价由 $x, y, v, \theta$与参考状态的误差项与$a, \omega$的二次项构成，终端代价由终端代价的误差项构成
- 约束
  - 等式约束：运动学模型（线性模型，简化对运动学模型的线性化处理）+ $\omega, v$到左右前轮转角的转换
  - 不等式约束：速度范围，加速度范围，左右前轮转角范围
$$
\begin{align*}
& \mathbf{J} = (x_N-x_{r,N})^TQ_N(x_N-x_{r,N}) + \sum_{i=0}^{N-1}{\{(x_i-x_{r,i})^TQ_i(x_i-x_{r,i})+u_i^TR_iu_i\}} \\
& \begin{align*}
    \mathrm{s.t.} &\\
    &\begin{align}
        & \qquad \begin{bmatrix} \dot{x} \\ \dot{y} \\ \dot{\theta} \\ \dot{v} \end{bmatrix} = \begin{bmatrix} v\cos{\theta} \\ v\sin{\theta} \\ \omega \\ a\end{bmatrix} \hspace{4cm} \\
        & \qquad \tan{\varphi_l} = \frac{2\omega l}{2v-b\omega} \\
        & \qquad \tan{\varphi_r} = \frac{2\omega l}{2v+b\omega} \\
        & \qquad l_{left} \le \varphi_l \le u_{left} \\
        & \qquad l_{right} \le \varphi_r \le u_{right} \\
        & \qquad v_{lb} \le v \le v_{ub} \\
        & \qquad a_{lb} \le a \le a_{ub} \\
        & \qquad \varphi_{lb} \le \varphi_{left} \le \varphi_{ub} \\
        & \qquad \varphi_{lb} \le \varphi_{right} \le \varphi_{ub} \\
        & \qquad \dot{\varphi}_{lb} \le \dot{\varphi_{left}} \le \dot{\varphi}_{ub} \\
        & \qquad \dot{\varphi}_{lb} \le \dot{\varphi_{right}} \le \dot{\varphi}_{ub} \\
    \end{align}
\end{align*}
\end{align*}
$$

### 1.3 线性化离散化
- 运动学约束
$$
\begin{align*}
& \begin{bmatrix} x_{k+1} \\ y_{k+1} \\ \theta_{k + 1} \\ v_{k + 1} \end{bmatrix} =
 (E+\Delta{t}\frac{\partial{f}}{\partial{x}}^T)x_k + \Delta{t}\frac{\partial{f}}{\partial{u}}^Tu_k +
 \Delta{t}(f(x_{r,k}, u_{r,k}) - \frac{\partial{f}}{\partial{x}}^Tx_{r,k}-\frac{\partial{f}}{\partial{u}}^Tu_{r,k}) \\
& \frac{\partial{f}}{\partial{x}} =
\begin{bmatrix}
0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 \\
-v_{r,k}\sin{\theta_{r,k}} & v_{r,k}\cos{\theta_{r,k}} & 0 & 0 \\
\cos{\theta_{t,k}} & \sin{\theta_{t,k}} & 0 & 0
\end{bmatrix} \\
& \frac{\partial{f}}{\partial{u}} =
\begin{bmatrix}
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1
\end{bmatrix}
\end{align*}
$$
- 前轮转角范围约束
$$
\begin{align*}
\begin{bmatrix}
0 & 0 & 0 & -2\tan{\varphi_{lb}} \\
0 & 0 & 0 & -2\tan{\varphi_{ub}} \\
0 & 0 & 0 & -2\tan{\varphi_{lb}} \\
0 & 0 & 0 & -2\tan{\varphi_{ub}} \\
\end{bmatrix} \mathrm{x_k} +
\begin{bmatrix}
-(2l-b\tan{\varphi_{lb}}) & 0 \\
(2l+b\tan{\varphi_{ub}}) & 0 \\
-(2l+b\tan{\varphi_{lb}}) & 0 \\
(2l-b\tan{\varphi_{ub}}) & 0 \\
\end{bmatrix}\mathrm{u_k} \le 0
\end{align*}
$$
- 前轮转角速度约束
采用前向差分近似计算转角速度，$\dot{\varphi}_{k} \approx \frac{\varphi_{k+1} - \varphi_{k}}{\Delta t}$
$$
\begin{align*}
& \tan{(\varphi_{k+1} - \varphi_{k})} = \frac{\tan\varphi_{k+1}-\tan\varphi_{k}}{1 + \tan\varphi_{k+1}\tan\varphi_{k}}
\end{align*}
$$
其中，$\tan\varphi$是$v,\omega$的函数，所以上式可以进一步表示为如下：
$$
\begin{align*}
& h_{k} = \tan{(\varphi_{k+1} - \varphi_{k})} = \frac{g_{k+1}-g_{k}}{1 + g_{k+1}g_{k}} \\
& \tan{(\Delta{t} \dot\varphi_{lb})} \le h_{k} \le \tan{(\Delta{t} \dot\varphi_{ub})}
\end{align*}
$$
将表达式在参考轨迹下泰勒展开，实现非线性约束的线性化
$$
p_k =
\begin{bmatrix}
\frac{\partial{h_{k}}}{\partial{v_{k}}} \\
\frac{\partial{h_{k}}}{\partial{\omega_{k}}} \\
\frac{\partial{h_{k}}}{\partial{v_{k+1}}} \\
\frac{\partial{h_{k}}}{\partial{\omega_{k+1}}} \\
\end{bmatrix}
= \begin{bmatrix}
-\frac{g_{k+1}^2\frac{\partial g_k}{\partial v_{k}} + \frac{\partial g_k}{v_k}}{(1 + g_{k+1}g_{k})^2} \\
-\frac{g_{k+1}^2\frac{\partial g_k}{\partial \omega_{k}} + \frac{\partial g_k}{\omega_k}}{(1 + g_{k+1}g_{k})^2} \\
\frac{g_{k}^2\frac{\partial g_{k+1}}{\partial v_{k+1}} + \frac{\partial g_{k+1}}{v_{k+1}}}{(1 + g_{k+1}g_{k})^2} \\
\frac{g_{k}^2\frac{\partial g_{k+1}}{\partial \omega_{k+1}} + \frac{\partial g_{k+1}}{\omega_{k+1}}}{(1 + g_{k+1}g_{k})^2} \\
\end{bmatrix}
$$
因此，得到如下线性化不等式约束：
$$
\begin{align*}
\tan{(\Delta{t} \dot\varphi_{lb})} \le
p_{k,r}^T \begin{bmatrix} v_k \\ \omega_k \\ v_{k + 1} \\ \omega_{k + 1}\end{bmatrix} + h_k^r - p_{k,r}^T\begin{bmatrix} v_{k,r} \\ \omega_{k,r} \\ v_{k + 1,r} \\ \omega_{k + 1,r}\end{bmatrix}
\le \tan{(\Delta{t} \dot\varphi_{ub})}
\end{align*}
$$