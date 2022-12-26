# Constrained-Optimization-Homework
Constrained optimization homework for the Numerical Optimization exam.

Text of the problem: consider the problem described by equation (3) in [1]. Use your own implementation of the projected gradient method to solve the problem with
$n = 10^d$ and $d = 3, 4, 5$, both using exact derivatives and finite differences to approximate the gradient. Compare the behavior of the two implementations, using the following values for the increment $h$:
```math
h = 10^{-k}\abs{\hat{x}}, \qquad k = 2, 4, 6, 8, 10, 12
```
where $\hat{x}$ is the point at which the derivatives have to be approximated.
