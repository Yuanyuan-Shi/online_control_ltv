# Stable Online Control of Linear Time-Varying Systems

This repository contains source code necessary to reproduce the results presented in the following paper: 
Stable Online Control of Linear Time-Varying Systems (http://proceedings.mlr.press/v144/qu21a/qu21a.pdf).

Authors: Guannan Qu, Yuanyuan Shi, Sahin Lale, Steven Low, Anima Anandkumar and Adam Wierman

Accepted and Presented at the 3rd Annual Conference on Learning for Dynamics and Control (L4DC).

### Paper Highlight

Linear time-varying (LTV) systems are widely used for modeling real-world dynamical systems
due to their generality and simplicity. Providing stability guarantees for LTV systems is one of the
central problems in control theory. However, existing approaches that guarantee stability typically
lead to significantly sub-optimal cumulative control cost in online settings where only current or
short-term system information is available. In this work, we propose an efficient online control
algorithm, COvariance Constrained Online Linear Quadratic (COCO-LQ) control, that guarantees
input-to-state stability for a large class of LTV systems while also minimizing the control cost. The
proposed method incorporates a state covariance constraint into the semi-definite programming
(SDP) formulation of the LQ optimal controller. We empirically demonstrate the performance of
COCO-LQ in both synthetic experiments and a power system frequency control example.

### Algorithm

![alt text](https://github.com/Yuanyuan-Shi/online_control_ltv/blob/main/figs/WX20220721-112402.png)


### Experiment

#### Synthetic Time-Varying Systems


#### Frequency Control with Renewable Generation

We consider a power system frequency control problem with time-varying renewable generation and thus, system inertia. The state space model 
of power system frequency dynamics follow, <img src="https://render.githubusercontent.com/render/math?">

$$
\begin{bmatrix}
\dot{\theta}\\
\dot{\omega} 
\end{bmatrix} = 
\begin{bmatrix}
0 & I \\
-M_t^{-1} L & -M_t^{-1} D
\end{bmatrix}
\begin{bmatrix}
\theta \\
\omega
\end{bmatrix} + 
\begin{bmatrix} 
0 \\ 
M_t^{-1}\end{bmatrix} P_{in}
$$

where the state variable is defined as the stacked vector of the voltage angle $\theta$ and frequency $\omega$.
$M_t = diag(m_{t,i})$ is the inertia matrix, where $m_{t,i}$ represents the equivalent rotational inertia at
bus i and time t. $M_t$ is time-varying and depends on the mix of online generators, since only
thermal generators provide rotational inertia and renewable generation does not. $D = diag(d_i)$ is the damping matrix, 
where $d_i$ is the generator damping coefficient. $L$ is the network susceptance matrix. The control variable pin 
corresponds to the electric power generation. 

We test on the  IEEE WECC 3-machine 9-bus system, where the system is changing between two states: a high renewable generation 
scenario where $m_{t,i} = 2$ (i.e., 80 percent renewable with zero inertia and 20 percent of thermal generation with 10s inertia), 
and a low renewable generation scenario where $m_{t,i} = 8$ (i.e., 20 percent renewable and 80 percent thermal generation), with 
additional random fluctuations between $[0, 0.2]$.

![alt text](https://github.com/Yuanyuan-Shi/online_control_ltv/blob/main/figs/WX20220721-112345.png)

