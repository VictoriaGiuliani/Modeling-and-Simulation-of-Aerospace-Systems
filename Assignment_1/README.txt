TASK DESCRIPTION ASSIGNMENT 1

- EXERCISE 1
Given the function f(x)= cos x - x, guess a and b such that f(a)f(b)<0. Find the zero(s) of f in [a,b] with 
     * the bisection method,
     * the secant method,
     * the regula falsi method.
 
All solutions must be calculated with 8-digit accuracy. Which method requires the least number of function evaluations? 
Report the computational time required by each method.

- EXERCISE 2
Let f be a two-dimensional vector-valued function f(x) = (x1^2 −x1 − x2; x1^2/16 + x2^2 − 1 )where x = (x1; x2).
Find the zero(s) of f by using Newton’s method with ∂f/∂x
1) computed analytically, and 
2) estimated through finite differences. Which version is more accurate?

- EXERCISE 3
The Initial Value Problem ̇x = x − t^2 + 1, x(0) = 12, has analytic solution x(t) = t^2 + 2t + 1 − 12e^t.
1) Implement a general-purpose, fixed-step Heun’s method (RK2);
2) Solve the IVP in t [0,2] for h1= 0.5, h2= 0.2, h3= 0.05, h4= 0.01 and compare the numerical vs the analytical solution; 
3) Repeat points 1)–2) with RK4;
4) Trade off between CPU time & integration error

- EXERCISE 4
Let ̇x = A(α)x be a two-dimensional system with A(α) = [0, 1; −1, 2cosα].  Notice that A(α) has a pair of complex conjugate eigenvalues
on the unit circle; α denotes the angle from the Re{λ}-axis.  
1) Write the operator FRK2(h,α) that maps x_k into x_k+1, namely x_k+1 = FRK2(h,α)x_k. 
2)With α = π, solve the problem “Find h ≥ 0 s.t. max (|eig(F(h,α))|) = 1”.
3) Repeat point 2) for α [0,π] and draw the solutions in the(hλ)-plane. 
4) Repeat points 1)–3) with RK4 and represent the points {hiλ} of Exercise 3 with t= 0. What can you say?

- EXERCISE 5
Consider the IVP ̇x = A(α)x, x(0) = [1,1]T, to be integrated in t∈[0,1]. 
1) Take α∈[0,2π] and solve the problem “Find h ≥ 0 s.t. ‖ xan(1)−xRK1(1) ‖∞= tol”, where xan(1) and xRK1(1) are the analytical 
and the numerical solution (with RK1) at the final time, respectively, and tol = {10−3, 10−4, 10−5, 10−6}. 
2) Plot the five locus of solutions in the (hλ)-plane; plot also the function evaluations vs tol for α = π. 
3) Repeat points 1)–2) for RK2 and RK4.

- EXERCISE 6
Consider the backinterpolation method BI20.4. 
1) Derive the expression of the linear operator BBI20.4(h,α) such that xk+1 = BBI20.4(h,α)xk. 
2) Following the approach of point 3) in Exercise5, draw the stability domain of BI20.4 in the(hλ)-plane. 
3) Derive the domain of numerical stability of BI2θ for the values of θ= [0.1,0.3,0.7,0.9].

- EXERCISE 7
Consider the IVP ̇x = Bx with B = [−180.5,219.5; 179.5,−220.5] and x(0) = [1,1]' to be integrated in t∈[0,5]. Notice that x(t) =e^(Bt) x(0). 
1) Solve the IVP using RK4 with h= 0.1;
2) Repeat point 1) using BI20.1; 
3) Compare the numerical results in points 1) and 2) against the analytic solution; 
4) Compute the eigenvalues associated to the IVP and represent them on the(hλ)-plane both for RK4 and BI20.1; 
5) Discuss the results.