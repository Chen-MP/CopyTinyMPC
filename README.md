 
# TinyMPC
[1] Nguyen K , Schoedel S , Alavilli A ,et al.TinyMPC: Model-Predictive Control on Resource-Constrained Microcontrollers[J].IEEE, 2023.DOI:10.1109/ICRA57147.2024.10610987.

--- 

## 1. Linear Quadratic Regulator(LQR)

### 1.1. Cost Function:
- 论文[1] 公式(1) 在有限时域内最小化如下二次代价函数：

$$
\begin{aligned}
\min_{x_{1:N}, \; u_{1:N-1}} J \; = &\; 
\frac{1}{2} x_N^\top Q_f x_N + q_f^\top x_N \\&
\;+\; \sum_{k=1}^{N-1} 
\Big( \frac{1}{2} x_k^\top Q x_k + q_k^\top x_k + \frac{1}{2} u_k^\top R u_k + r_k^\top u_k \Big) \tag{1.1}
\end{aligned}
$$

- 系统模型等式约束 && 状态变量、输入变量的不等式约束：

    $$ 
 \text{subject to} \quad x_{k+1} = A x_k + B u_k, \quad \forall k \in \{1, 2, \dots, N-1\}, \\
 x_k \in \mathcal{X}, \quad u_k \in \mathcal{U}.
$$
    下面写出上式的推导过程。  

### 1.2. 逐项求和形式：

- 预测域长度 $N$，状态 $x_k$ {$k=1,\dots,N$}，控制 $u_k${$k=1,\dots,N-1$}，参考轨迹 $x_k^{\mathrm{ref}},u_k^{\mathrm{ref}}$。
Cost Function 原式如下： 

    $$
\begin{aligned}
J \; = &\; \frac{1}{2} (x_N - x_N^{\mathrm{ref}})^\top Q_f (x_N - x_N^{\mathrm{ref}}) \\& 
\; + \; \frac{1}{2} \sum_{k = 1}^{N-1} \Big[ (x_k - x_k^{\mathrm{ref}})^\top Q (x_k - x_k^{\mathrm{ref}})+ (u_k - u_k^{\mathrm{ref}})^\top R (u_k - u_k^{\mathrm{ref}}) \Big] \tag{1.2}
\end{aligned}
$$

    其中 $Q \succeq 0,\; Q_f \succeq 0,\; R \succ 0$。

- 展开后（省略常数项）：

    $$
\begin{aligned}
J \; = &\; \frac{1}{2} x_N^\top Q_f x_N - (Q_f x_N^{\mathrm{ref}})^\top x_N \\&
\;+\; \sum_{k=0}^{N-1} \Big[ \frac{1}{2} x_k^\top Q x_k - (Q x_k^{\mathrm{ref}})^\top x_k+ \frac{1}{2} u_k^\top R u_k - (R u_k^{\mathrm{ref}})^\top u_k \Big] \;+\; \text{const.} 
\tag{1.3}
\end{aligned}
$$

- 将 式(1.1) 与 式(1.3) 对比可得：
式(1.1) 中一阶线性项的表达式
$$
\mathbf{q} =
\begin{bmatrix}
q_1 \\ q_2 \\ \vdots \\ q_{N-1} \\ q_f
\end{bmatrix}, 
\quad
\mathbf{r} =
\begin{bmatrix}
r_1 \\ r_2 \\ \vdots \\ r_{N-1}
\end{bmatrix}
$$

    并且通过参考轨迹得到：

$$
q_k = - Q x_k^{\mathrm{ref}}, \quad r_k = - R u_k^{\mathrm{ref}}, \quad q_f = - Q_f x_N^{\mathrm{ref}}
$$

- 此外补充说明：

    **① 状态变量$ \mathbf{x} $ && 控制变量$ \mathbf{u} $ 矩阵形式**
    $$
\mathbf{x} =
\begin{bmatrix}
x_1 \\ x_2 \\ \vdots \\ x_N
\end{bmatrix} \in \mathbb{R}^{N n_x}, 
\quad
\mathbf{u} =
\begin{bmatrix}
u_1 \\ u_2 \\ \vdots \\ u_{N-1}
\end{bmatrix} \in \mathbb{R}^{(N-1) n_u}
$$

    **② 参考轨迹矩阵形式：**

    $$
\mathbf{x}^{\mathrm{ref}} =
\begin{bmatrix}
x_1^{\mathrm{ref}} \\ x_2^{\mathrm{ref}} \\ \vdots \\ x_N^{\mathrm{ref}}
\end{bmatrix}, \quad
\mathbf{u}^{\mathrm{ref}} =
\begin{bmatrix}
u_1^{\mathrm{ref}} \\ \vdots \\ u_{N-1}^{\mathrm{ref}}
\end{bmatrix}
$$

    **③ 权重对角矩阵：**

    $$
\mathbf{Q} = \mathrm{diag}(Q, Q, \dots, Q, Q_f) \in \mathbb{R}^{N n_x \times N n_x}
$$

    $$
\mathbf{R} = \mathrm{diag}(R, R, \dots, R) \in \mathbb{R}^{(N-1) n_u \times (N-1) n_u}
$$


### 1.3. Riccati Recursion 黎卡提递推

- **终端条件：**
    $$
V_N(x_N) = \tfrac{1}{2} x_N^\top Q_f x_N + q_f^\top x_N
$$

    令

    $$
P_N = Q_f, \quad p_N = q_f
$$

- ***第一步：写出一般形式的 Value Function***
    $$
V_k(x_k) = \tfrac{1}{2} x_k^\top P_k x_k + p_k^\top x_k + c_k
$$

- ***第二步：在 $k = N-1$ 时展开***
    $$
    \begin{aligned}
V_{N-1}(x_{N-1}) \; =& \; \tfrac{1}{2} x_{N-1}^\top Q x_{N-1} + q_{N-1}^\top x_{N-1}+ \tfrac{1}{2} u_{N-1}^\top R u_{N-1} + r_{N-1}^\top u_{N-1} \\& + V_N(x_N)
\end{aligned}
$$
    又因为：
    $$
x_N = A x_{N-1} + B u_{N-1}
$$

    将上式代入 Value Function $ V_{N-1}(x_{N-1}) $ 可得：

    $$
    \begin{aligned}
V_{N-1}(x_{N-1}) \; =& \; \tfrac{1}{2} x_{N-1}^\top Q x_{N-1} + q_{N-1}^\top x_{N-1}+ \tfrac{1}{2} u_{N-1}^\top R u_{N-1} + r_{N-1}^\top u_{N-1} \\& + \tfrac{1}{2} (A x_{N-1} + B u_{N-1})^\top P_N (A x_{N-1} + B u_{N-1}) + p_N^\top (A x_{N-1} + B u_{N-1})
\end{aligned}
$$


- **第三步：对 $u_{N-1}$ 求导并令其为零**
    $$
\frac{\partial V_{N-1}}{\partial u_{N-1}} =
R u_{N-1} + r_{N-1} + B^\top P_N (A x_{N-1} + B u_{N-1}) + B^\top p_N = 0
$$

- **第四步：解出最优控制律**
    $$
(R + B^\top P_N B) u_{N-1} =- B^\top P_N A x_{N-1} - B^\top p_N - r_{N-1}
$$

    $$
u_{N-1}^\star = -K_{N-1} x_{N-1} - d_{N-1}
$$

    其中：

    $$
K_{N-1} = (R + B^\top P_N B)^{-1} B^\top P_N A
$$

    $$
d_{N-1} = (R + B^\top P_N B)^{-1} (B^\top p_N + r_{N-1})
$$

- **第五步：得到新的 Value Function 系数**
    **① 首先将最优控制律 $u^* $ 代回 Value Function $V_{N-1}$**

    $$
    \begin{aligned}
V_{N-1}(x_{N-1}, u_{N-1}) \; =& \; \tfrac{1}{2} x_{N-1}^\top Q x_{N-1} + q_{N-1}^\top x_{N-1} + \tfrac{1}{2} u_{N-1}^\top R u_{N-1} + r_{N-1}^\top u_{N-1} + \\& 
 V_N(A x_{N-1} + B u_{N-1})
\end{aligned}
$$

    得：

    $$
\begin{aligned}
V_{N-1}(x_{N-1}) &= \tfrac{1}{2} x_{N-1}^\top Q x_{N-1} + q_{N-1}^\top x_{N-1} \\
&\quad + \tfrac{1}{2} (-K_{N-1} x_{N-1} - d_{N-1})^\top R (-K_{N-1} x_{N-1} - d_{N-1}) \\
&\quad + r_{N-1}^\top (-K_{N-1} x_{N-1} - d_{N-1}) \\
&\quad + V_N\big(A x_{N-1} + B(-K_{N-1} x_{N-1} - d_{N-1})\big)
\end{aligned}
$$

    **② 展开每一项**
    控制代价展开：
    $$
\tfrac{1}{2} u^\top R u = \tfrac{1}{2} x_{N-1}^\top K_{N-1}^\top R K_{N-1} x_{N-1} + x_{N-1}^\top K_{N-1}^\top R d_{N-1} + \tfrac{1}{2} d_{N-1}^\top R d_{N-1}
$$

    线性项：
    $$
r_{N-1}^\top u = - x_{N-1}^\top K_{N-1}^\top r_{N-1} - r_{N-1}^\top d_{N-1}
$$

    递推项：
    $$
    \begin{aligned}
V_N(A x_{N-1} + B u) \; =& \; \tfrac{1}{2} \Big( (A - B K_{N-1}) x_{N-1} + B d_{N-1} \Big) ^\top P_N \Big( (A - B K_{N-1}) x_{N-1} + B d_{N-1} \Big) \\& + p_N^\top \Big( (A - B K_{N-1}) x_{N-1} + B d_{N-1} \Big)
\end{aligned}
$$

    **③ 整理成标准二次型**
    $$
V_{N-1}(x_{N-1}) = \tfrac{1}{2} x_{N-1}^\top P_{N-1} x_{N-1} + p_{N-1}^\top x_{N-1} + \text{const}
$$

    **④ 可得系数为：**

    $$
\begin{aligned}
P_{N-1} &= Q + K_{N-1}^\top R K_{N-1} + (A - B K_{N-1})^\top P_N (A - B K_{N-1}) \\
p_{N-1} &= q_{N-1} + (A - B K_{N-1})^\top (p_N - P_N B d_{N-1}) + K_{N-1}^\top (R d_{N-1} - r_{N-1})
\end{aligned}
$$

- **第六步：推广到一般 $k$**
对于任意 $k = N-1, N-2, \dots, 1$，递推得到：

    $$
u_k^* = - K_k  x_k - d_k
$$

    $$
K_k = (R + B^\top P_{k+1} B)^{-1} B^\top P_{k+1} A
$$

    $$
d_k = (R + B^\top P_{k+1} B)^{-1} (B^\top p_{k+1} + r_k)
$$

    $$
P_k = Q + K_k^\top R K_k + (A - B K_k)^\top P_{k+1} (A - B K_k)
$$

    $$
p_k = q_k + (A - B K_k)^\top (p_{k+1} - P_{k+1} B d_k) + K_k^\top (R d_k - r_k)
$$

**上面的式子与 论文[1] 式(3) 完全相同。**

---

## 2. The Alternating Direction Method of Multipliers (ADMM)

### 2.1. 重构优化问题
- **① 假设存在原始优化问题：**

    $$
\min_{x} \;\; f(x) \quad \text{s.t. } \; x \in \mathcal{C}
$$

    目标是找到 $ x^* $ 使得 $ f(x) $ 最小 ，且 $ x^* $ 必须属于凸集 $ \mathcal{C} $

- **② 引入松弛变量(Slack Variable) $ z $ && 指示函数(Indicator Function)**

    $$
I_{\mathcal{C}}(z) = 
\begin{cases}
0, & \text{if} \; z \in \mathcal{C} \\
+\infty, & \text{otherwise}
\end{cases}
$$

    于是原始优化问题的不等式约束问题转化为无约束问题：
$$
\min_{x, \; z} \;\; f(x) + I_{\mathcal{C}}(z) \quad \text{subject to} \;\; x = z
$$

- **③ 原理**
    **Queation:** 为什么等价？
    **Answer:** 
    原因是 约束 $ x = z $ 强制 $ x $ 只能取 $ z $ 的值，而 $ z $ 受 $ I_{\mathcal{C}}(z) $  约束：
    i. 当 $ z \in \mathcal{C} $ 时， $ I_{\mathcal{C}}(z) = 0 $，问题就变成了 $ \min_x f(x) $ 但要求 $ x=z $ ，也就是 $ x \in \mathcal{C} $。
    ii. 当 $ z \notin \mathcal{C} $ 时， $ I_{\mathcal{C}}(z) = \infty $，这会使得目标函数变成无穷大，因此优化器不会选取这样的 $ z $，即它隐式地强制了 $ z \in \mathcal{C} $。
    iii. 所以，优化器只能选择满足 $ x=z $ 且 $ z \in \mathcal{C} $ 的解，这就等价于原来的约束 $ x \in \mathcal{C} $。   
    **Queation:** 为什么要这样做？
    **Answer:**
    i. 原问题直接带有约束 $ x \in \mathcal{C} $，而在很多优化方法（如拉格朗日对偶、投影梯度法、ADMM）中，处理不等式约束比较麻烦。
    通过引入指示函数 $ I_{\mathcal{C}}(z) $ 和松弛变量 $ z $ ，不等式约束变成了等式约束 $ x=z $，这样就可以使用拉格朗日方法或增广拉格朗日法来求解。
    ii. 此外，不等式约束 $ x \in \mathcal{C} $ 可能难以直接优化，而等式约束 $ x=z $ 可以转换成投影运算：$ z^{(k+1)} = \operatorname{Proj}_{\mathcal{C}}(x^{(k)}) $，这样可以使用投影梯度法或交替方向乘子法（ADMM）。

### 2.2. Augmented Lagrange Function 增广拉格朗日函数

在 2.1 中将不等式约束问题转换为了等式约束问题，现在开始使用拉格朗日将等式约束问题转化为无约束问题。

- **1. 标准拉格朗日（无约束鞍点形式）**
    首先从以下等式约束问题出发：
    $$
\min_{x,z}\; f(x)+I_{\mathcal C}(z)\quad \text{subject to } x=z
$$

    引入**拉格朗日乘子**$ \lambda $并且构造拉格朗日函数：
    $$
\mathcal L(x,z,\lambda)=f(x)+I_{\mathcal C}(z)+\lambda^\top(x-z).
$$

    于是原问题等价于无约束问题：
    $$
\min_{x,z}\ \max_{\lambda}\ \mathcal L(x,z,\lambda).
$$

    这里 $ \lambda $ 是拉格朗日乘子，作用于等式约束 $ x=z $ 。
    好处在于将等式约束合并到目标函数，转化为无约束优化问题。
    
- **2. 增广拉格朗日**

    为了改善数值性质，减少数值震荡，引入二次罚项，得到增广拉格朗日：
    $$
\mathcal L_\rho(x,z,\lambda) = f(x)+I_{\mathcal C}(z) + \lambda^\top(x-z)
+\frac{\rho}{2}\|x-z\|^2, 
\qquad \rho>0.
$$

    **此时，上式完全与 论文[1] 公式(9) 完全相同。**

    同样是无约束鞍点问题：
    $$
\min_{x,z}\ \max_{\lambda}\ \mathcal L_\rho(x,z,\lambda).
$$

### 2.3. ADMM

- **1. ADMM 更新步骤与理解**

    在第 $k$ 次迭代中，ADMM 包含以下三个更新步骤：

    **① primal update $x$：** 固定 $z, \lambda$，单独优化 $x$
    
    $$
  x^{k+1} = \arg\min_x \; \mathcal{L}_\rho(x, z^k, \lambda^k)
  $$

    **② slack update $z$：** 固定 $x, \lambda$，单独优化 $z$
    $$
  z^{k+1} = \arg\min_z \; \mathcal{L}_\rho(x^{k+1}, z, \lambda^k)
  $$

    **③ dual update $\lambda$：** 通过更新对偶变量，逐渐强化对约束的惩罚力度
    $$
  \lambda^{k+1} = \lambda^k + \rho \big(x^{k+1} - z^{k+1})
  $$

- **2. 关键直觉**：  
  $x$ 和 $z$ 分别针对自己部分的问题优化；  
  $\lambda$ 确保约束逐步被满足；  
  这样通过交替更新，原本耦合的难题被拆分成了多个容易的小问题。

---

## 3. Combining LQR and ADMM

**综合 第一部分 LQR 和 第二部分 ADMM 可得：**

- **1. 定义松弛变量：** $x_k=z_k,\; u_k=w_k$。

- **2. 定义指示函数：** $I_{\mathcal X}(z)$; $I_{\mathcal U}(w)$
  
- **3. 重构代价函数得到新的优化问题：**

    $$ 
\begin{aligned} 
\mathcal L \; = &\; \frac{1}{2} (x_N - x_N^{\mathrm{ref}})^\top Q_f (x_N - x_N^{\mathrm{ref}}) \\& 
\; + \; \frac{1}{2} \sum_{k = 1}^{N-1} \Big[ (x_k - x_k^{\mathrm{ref}})^\top Q (x_k - x_k^{\mathrm{ref}})+ (u_k - u_k^{\mathrm{ref}})^\top R (u_k - u_k^{\mathrm{ref}}) \Big] \\&
\; + \; I_{\mathcal X}(z) \; + \; I_{\mathcal U}(w) \\&
\; + \; \sum_{k=1}^{N} \Big[ \lambda_k^\top (x_k - z_k) \; + \; \frac{\rho}{2} (x_k - z_k)^\top (x_k - z_k) \Big] \\&
\; + \; \sum_{k=1}^{N-1} \Big[ \mu_k^\top (u_k - w_k) \; + \; \frac{\rho}{2} (u_k - w_k)^\top (u_k - w_k) \Big] \\&
\text{subject to } \; x_{k+1} = A x_{k} + B u_{k}
 \end{aligned} 
 $$

### 3.1. 使用 ADMM 三步更新：

$$
\begin{aligned}
\min_{x,u} \quad & \frac{1}{2}(x_N - x_N^{\mathrm{ref}})^\top Q_f (x_N - x_N^{\mathrm{ref}}) \\
& + \frac{1}{2}\sum_{k=1}^{N-1} \left[ (x_k - x_k^{\mathrm{ref}})^\top Q (x_k - x_k^{\mathrm{ref}}) + (u_k - u_k^{\mathrm{ref}})^\top R (u_k - u_k^{\mathrm{ref}}) \right] \\
& + \sum_{k=1}^{N} \left[ \lambda_k^\top (x_k - z_k) + \frac{\rho}{2} \|x_k - z_k\|_2^2 \right] \\
& + \sum_{k=1}^{N-1} \left[ \mu_k^\top (u_k - w_k) + \frac{\rho}{2} \|u_k - w_k\|_2^2 \right] \\
& +  I_{\mathcal X}(z) \; + \; I_{\mathcal U}(w) \\
\text{s.t.} \quad & x_{k+1} = A x_k + B u_k, \quad k=1,\dots,N
\end{aligned}
$$

#### 3.1.1. 第一步：更新主要变量 $(x, u)$

- **① 将上面的增广拉格朗日函数重新整理：**

    对于状态 $x_k$（$k=1,\dots,N-1$）：
    
    $$
\frac{1}{2}x_k^\top Q x_k - x_k^\top Q x_k^{\mathrm{ref}} + \lambda_k^\top x_k + \frac{\rho}{2}x_k^\top x_k - \rho x_k^\top z_k
$$

    合并：
    
    $$
\frac{1}{2}x_k^\top (Q + \rho I) x_k + (-Q x_k^{\mathrm{ref}} + \lambda_k - \rho z_k)^\top x_k
$$

    对于终端状态 $x_N$：
    
    $$
\frac{1}{2}x_N^\top Q_f x_N - x_N^\top Q_f x_N^{\mathrm{ref}} + \lambda_N^\top x_N + \frac{\rho}{2}x_N^\top x_N - \rho x_N^\top z_N
$$

    合并：

    $$
\frac{1}{2}x_N^\top (Q_f + \rho I) x_N + (-Q_f x_N^{\mathrm{ref}} + \lambda_N - \rho z_N)^\top x_N
$$

    对于控制输入 $u_k$（$k=1,\dots,N-1$）：
    
    $$
\frac{1}{2}u_k^\top R u_k - u_k^\top R u_k^{\mathrm{ref}} + \mu_k^\top u_k + \frac{\rho}{2}u_k^\top u_k - \rho u_k^\top w_k
$$

    合并：
    
    $$
\frac{1}{2}u_k^\top (R + \rho I) u_k + (-R u_k^{\mathrm{ref}} + \mu_k - \rho w_k)^\top u_k
$$

    **综上可得，定义修改的权重矩阵和向量：**

    $$
 \begin{aligned}
     \tilde{Q}_k &= Q + \rho I \\
     \tilde{Q}_f &= Q_f + \rho I \\
     \tilde{R}_k &= R + \rho I \\
     \tilde{q}_k &= -Q x_k^{\mathrm{ref}} + \lambda_k - \rho z_k = q_k + \lambda_k - \rho z_k \\
     \tilde{q}_f &= -Q_f x_N^{\mathrm{ref}} + \lambda_N - \rho z_N = q_f + + \lambda_N - \rho z_N \\
     \tilde{r}_k &= -R u_k^{\mathrm{ref}} + \mu_k - \rho w_k = r_k + \mu_k - \rho w_k
     \end{aligned}
 $$

    **与 论文[1] 公式(15) 相同。重新整理目标函数可得：**

    $$
\begin{aligned}
\min_{x,\,u} \quad 
& \tfrac{1}{2} x_N^\top \tilde{Q}_f x_N + \tilde{q}_f^\top x_N + \sum_{k=1}^{N-1} \Big[ \tfrac{1}{2} x_k^\top \tilde{Q}_k x_k + \tilde{q}_k^\top x_k + \tfrac{1}{2} u_k^\top \tilde{R}_k u_k + \tilde{r}_k^\top u_k \Big] \$$6pt]
\text{s.t.}\quad 
& x_{k+1} = A x_k + B u_k, \quad k=1,\dots,N
\end{aligned} \tag{3.1}
$$

    其中参数定义为：

    $$
\begin{aligned}
\tilde{Q}_k &= Q + \rho I, \quad & q_k &= -Q x_k^{\mathrm{ref}}, 
& \tilde{q}_k &= q_k + \lambda_k - \rho z_k, \$$6pt]
\tilde{Q}_f &= Q_f + \rho I, \quad & q_f &= -Q_f x_N^{\mathrm{ref}}, 
& \tilde{q}_f &= q_f + \lambda_N - \rho z_N, \$$6pt]
\tilde{R}_k &= R + \rho I, \quad & r_k &= -R u_k^{\mathrm{ref}}, 
& \tilde{r}_k &= r_k + \mu_k - \rho w_k.
\end{aligned}
$$

- **② 引入缩放对偶变量**

    定义缩放后的对偶变量为：

    $$
y_k = \tfrac{1}{\rho}\lambda_k, \qquad g_k = \tfrac{1}{\rho}\mu_k.
$$

    因此移项可得：

    $$
\lambda_k = \rho y_k, \qquad \mu_k = \rho g_k.
$$

- **③ 重写线性项**

    代入缩放对偶变量后，可以将修正的线性项写为更紧凑的形式：

    $$
\begin{aligned}
\tilde{q}_k &= q_k + \lambda_k - \rho z_k 
= q_k + \rho y_k - \rho z_k 
= q_k + \rho (y_k - z_k), \$$6pt]
\tilde{q}_f &= q_f + \rho (y_N - z_N), \$$6pt]
\tilde{r}_k &= r_k + \rho (g_k - w_k).
\end{aligned}
$$

- **④ 黎卡提递推，更新主要变量**

    **观察上面的式子 (3.1)，可知已经重新转化为标准的有限时域 LQR 问题，使用动态规划求解。**

    **i. Backpropagation (Calculating Gain Matrix) 反向传播 计算增益矩阵**

    $$
K_k = (\tilde{R} + B^\top P_{k+1} B)^{-1} B^\top P_{k+1} A
$$

    $$
d_k = (\tilde{R} + B^\top P_{k+1} B)^{-1} (B^\top p_{k+1} + \tilde{r}_k)
$$

    $$
P_k = \tilde{Q} + K_k^\top \tilde{R} K_k + (A - B K_k)^\top P_{k+1} (A - B K_k)
$$

    $$
p_k = \tilde{q}_k + (A - B K_k)^\top ( \tilde{p}_{k+1} - P_{k+1} B d_k) + K_k^\top (R d_k - \tilde{r}_k)
$$

    **ii. Forward Deduction (Calculating Optimal Control) 前向传播 计算最优控制**

    $$
u_k^* = - K_k  x_k^* - d_k
$$

    $$
x_{k+1}^* = A x_k^* + B u_k^*
$$

    **下面 $t$ 表示迭代次数**
    **得到 预测控制序列  \{ $u_1^{t+1}, u_2^{t+1}, \dots, u_{N-1}^{t+1}$ \}**
    **和  预测状态轨迹 \{ $x_1^{t+1}, x_2^{t+1}, \dots, x_{N-1}^{t+1}$ \}**

#### 3.1.2. 更新松弛变量 $(z, w)$

固定 $x, u, \lambda, \mu$，最小化：
$$
\min_{z,w} \quad I_{\mathcal{X}}(z) + I_{\mathcal{U}}(w) + \sum_{k=1}^{N} \left[ \lambda_k^\top (x_k - z_k) + \frac{\rho}{2} \|x_k - z_k\|_2^2 \right] + \sum_{k=1}^{N-1} \left[ \mu_k^\top (u_k - w_k) + \frac{\rho}{2} \|u_k - w_k\|_2^2 \right]
$$

**分解为独立问题：**

对于每个 $z_k$：
$$
z_k^{t+1} = \arg\min_{z_k \in \mathcal{X}} \left[ \lambda_k^\top (x_k - z_k) + \frac{\rho}{2} \|x_k - z_k\|_2^2 \right] = \Pi_{\mathcal{X}} \left( x_k^{t+1} + y_k^t \right)
$$

对于每个 $w_k$：
$$
w_k^{t+1} = \arg\min_{w_k \in \mathcal{U}} \left[ \mu_k^\top (u_k - w_k) + \frac{\rho}{2} \|u_k - w_k\|_2^2 \right] = \Pi_{\mathcal{U}} \left( u_k^{t+1} + g_k^t \right)
$$

#### 3.1.3. 更新乘子 对偶变量 $(\lambda, \mu)$

乘子通过梯度上升更新：
$$
\lambda_k^{t+1} = \lambda_k + \rho (x_k^{t+1} - z_k^{t+1})
$$

$$
\mu_k^{t+1} = \mu_k + \rho (u_k^{t+1} - w_k^{t+1})
$$

又  $ \quad \lambda_k = \rho y_k, \quad \mu_k = \rho g_k $ 可得：

$$
y_k^{t+1} = y_k^t + x_k^{t+1} - z_k^{t+1}
$$

$$
g_k^{t+1} = g_k^t + u_k^{t+1} - w_k^{t+1}
$$

---

## 4. Code Practice

### 4.1. ADMM 第一步 Primal Update

$$
K_k = (\tilde{R} + B^\top P_{k+1} B)^{-1} B^\top P_{k+1} A
$$

$$
d_k = (\tilde{R} + B^\top P_{k+1} B)^{-1} (B^\top p_{k+1} + \tilde{r}_k)
$$

$$
P_k = \tilde{Q} + K_k^\top \tilde{R} K_k + (A - B K_k)^\top P_{k+1} (A - B K_k)
$$

$$
p_k = \tilde{q}_k + (A - B K_k)^\top ( \tilde{p}_{k+1} - P_{k+1} B d_k) + K_k^\top (R d_k - \tilde{r}_k)
$$

**1. 首先是在动态规划过程中迭代 Value Function -> Hessian Matrix $P_k$ && Control Law -> Feedback Gain Matrix $K_k$，得到 $ P_{\mathrm{inf}} \; and \; K_{\mathrm{inf}} $**

**2. 然后发现可得 Precomputed Constant Matrices 预计算常量矩阵，计算缓存**

$$
C_1 = (\tilde{R} + B^\top P_{\mathrm{inf}} B)^{-1}
$$

$$
C_2 = (A − B K_{\mathrm{inf}} )^\top
$$

**3. 计算 预测时域（Prediction Horizon）中的 Value Function -> Gradients $p_k$ && Control Law -> Feedforward Vector (Compensation Term) $d_k$**
**4. 综合上面的步骤将会得到 预测状态轨迹 \{ $ x_1, x_2, \; \dots \; , x_N $ \} and 预测控制序列 \{ $ u_1, u_2, \; \dots \; , u_{N-1} $ \}**

### 4.2. ADMM 第二步 Slack Update
### 4.3. ADMM 第三步 Dual Update
上述两步都是与公式相同，不再赘述。

### 4.4. 对 $q_k$ 简化的推导过程

- **首先我们有理论公式：**

    $$
p_k = q_k + (A - B K_k)^\top (p_{k+1} - P_{k+1} B d_k) + K_k^\top (R d_k - r_k)
$$


    $$
p_k = q_k + (A - B K_k)^\top p_{k+1} - (A - B K_k)^\top P_{k+1} B d_k + K_k^\top R d_k - K_k^\top r_k
$$

    需要证明：

    $$
-(A - B K_k)^\top P_{k+1} B d_k + K_k^\top R d_k = 0
$$

- **第一步：展开括号**

    $$
-(A - B K_k)^\top P_{k+1} B d_k = -A^\top P_{k+1} B d_k + K_k^\top B^\top P_{k+1} B d_k
$$

    所以原式变为：

    $$
-A^\top P_{k+1} B d_k + K_k^\top B^\top P_{k+1} B d_k + K_k^\top R d_k
$$

  提取公因子：

    $$
\text{上式} = -A^\top P_{k+1} B d_k + K_k^\top (B^\top P_{k+1} B + R) d_k
$$

- **第二步：代入 $K_k$ 的定义**

    根据 Riccati 方程，我们知道：

    $$
K_k = (R + B^\top P_{k+1} B)^{-1} B^\top P_{k+1} A = Q_{uu}^{-1} B^\top P_{k+1} A
$$

    所以：

    $$
K_k^\top = A^\top P_{k+1} B Q_{uu}^{-1}
$$

- **第三步：代入 $K_k^\top$**

    $$
\begin{aligned}
& -A^\top P_{k+1} B d_k + K_k^\top (B^\top P_{k+1} B + R) d_k \\
=& -A^\top P_{k+1} B d_k + (A^\top P_{k+1} B Q_{uu}^{-1}) (Q_{uu}) d_k \\
=& (-A^\top P_{k+1} B d_k) + (A^\top P_{k+1} B d_k) \\
=& 0
\end{aligned}
$$

    **完美抵消**

- **完整的推导结果**

    因此：

    $$
\begin{aligned}
p_k =& q_k + (A - B K_k)^\top p_{k+1} - \Big[ (A - B K_k)^\top P_{k+1} B d_k + K_k^\top R d_k \Big] - K_k^\top r_k \\
=& q_k + (A - B K_k)^\top p_{k+1} - \Big[0 \Big] - K_k^\top r_k \\
=& q_k + (A - B K_k)^\top p_{k+1} - K_k^\top r_k
\end{aligned}
$$

最终公式化简为：

$$
p_k 
$$

这正是代码中的实现形式！
