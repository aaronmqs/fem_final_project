# fem_final_project
Repository with all shared codes for the final project of Finite Element Methods.
The notation used in the code follows what's described in the [Project](https://github.com/aaronmqs/fem_final_project/blob/main/project.pdf) and in the [notation](https://github.com/aaronmqs/fem_final_project/blob/main/notation.m).

# Questions:

- [ ] Hollow cylinder: the cylinder gets too weird when applying a force.

# What was done so far:

- [x] Part 1
  - [x] Part 1.1
    - [x] 1
    - [x] 2
  - [x]  Part 1.2
    - [x] 1
    - [x] 2
    - [x] 3
    - [x] 4
  - [x] Part 1.3
    - [x] 1
    - [x] 2
    - [x] 3
    - [x] 4
    - [x] 5
    - [x] 6
  - [x] Part 1.4
    - [x] 1
    - [x] 2
    - [x] 3
    - [x] 4 [DUE TO 04/05] (changed to 04/08 in class)
- [x] Part 2
  - [x] Part 2.1
  - [x] Part 2.2
    - [x] 1
    - [x] 2
  - [x] Part 2.3
    - [x] 1
    - [x] 2
  - [x] Part 2.4 
    - [x] 1
    - [x] 2
- [x] Part 3
  - [x] Part 3.1
  - [x] Part 3.2 [DUE TO 04/19]
- [x] Part 4
    - [x] 1
    - [x] 2
    - [x] 3
    - [x] 4
    - [x] 5
    - [x] 6
- [x] Part 5
  - [x] 1
  - [x] 2
  - [x] 3
  - [x] 4
- [ ] Part 6
  - [x] 1
  - [x] 2
- [x] Part 7 
  - [x] 1
  - [x] 2 [DUE TO 05/01]

# Answered questions

- [x] The number of quadrature nodes for simplices: see Task 5 for Part 1.3.
      
      -- check_moments_refdom("simp", 2, 1)
      Whe running this function, I need more than ceil((porder + 1) / 2) quadrature points per dimension in order to get accurate moments for the simplex. Why?

      --Answer:
      This happens because the simplex quadrature nodes are obtained via a transformation from the hypercube quadrature nodes. This transformation "increases the polynomial degree" of the integrand, meaning that more quadrature nodes would be required for the exact integration.
      
- [x] See Task 6 for Part 1.3 (simplex).
      
      -- driver
      The errors don't decrease even with the increase of quadrature points. Is it because the added qpoints are accumulating in one place?

      -- Answer:
      The error actually decreases: the problem is that the "exact value" I'm using to compare was provided with "only" four decimal digits, meaning that it's not the exact value of the integral.

- [x] See Task 6 for Part 1.3 (hcube).

      -- driver
      Same issue: more qpoints than expected.

      -- Answer:
      Even if the integrals are evaluated at a hypercube element, this is not an idealized element. This means that the integrand is actually a composition of functions from the idealized domain to the physical domain. This transformation "increases the polynomial degree" (quoted because the transformation and other new terms included in the integrand can even not be a polynomial), which means that more quadrature points are needed.

- [x] Simplicial mesh from hcube mesh (see Task 2 for Part 1.4.)

      -- The same function that creates hypercube element meshes also creates the simplicial meshes.
      
- [x] Part 4.6: why is the optimal convergence rate for poder = 2 so high? Is it because there's a 2nd degree polynomial as part of the exact solution?

      -- It's ok to have an optimal convergence rate higher than expected. My ituition is that the solution is composed by a quadratic polynomial, and therefore p=2 is kind of good in this case.

- [x] Part 5.3: Is [~, f2v, ~] = create_nodes_bndy_refdom_simp(ndim, porder); the same as msh.lfcnsp.f2v? (Same thing on Part 5.4)

      -- Yeah, they're exactly the same.
