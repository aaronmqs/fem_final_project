# fem_final_project
Repository with all shared codes for the final project of Finite Element Methods.
The notation used in the code follows what's described in the [Project](https://github.com/aaronmqs/fem_final_project/blob/main/project.pdf) and in the [notation](https://github.com/aaronmqs/fem_final_project/blob/main/notation.m).

# Questions

- [ ] The number of quadrature nodes for simplices: see Task 5 for Part 1.3.
      
      -- check_moments_refdom("simp", 2, 1)
      Whe running this function, I need more than ceil((porder + 1) / 2) quadrature points per dimension in order to get accurate moments for the simplex. Why?
      
- [ ] See Task 6 for Part 1.3 (simplex).
      
      -- driver
      The errors don't decrease even with the increase of quadrature points. Is it because the added qpoints are accumulating in one place?

- [ ] See Task 6 for Part 1.3 (hcube).

      -- driver
      Same issue: more qpoints than expected.

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
  - [ ] Part 1.4
    - [x] 1
    - [ ] 2
    - [ ] 3
    - [ ] 4 [DUE TO 04/05] (changed to 04/08 in class)
- [ ] Part 2
  - [ ] Part 2.1
  - [ ] Part 2.2
    - [ ] 1
    - [ ] 2
  - [ ] Part 2.3
    - [ ] 1
    - [ ] 2
  - [ ] Part 2.4 
    - [ ] 1
    - [ ] 2
- [ ] Part 3
  - [ ] Part 3.1
  - [ ] Part 3.2 [DUE TO 04/19]
- [ ] Part 4
    - [ ] 1
    - [ ] 2
    - [ ] 3
    - [ ] 4
    - [ ] 5
    - [ ] 6
- [ ] Part 5
  - [ ] 1
  - [ ] 2
  - [ ] 3
  - [ ] 4
- [ ] Part 6
  - [ ] 1
  - [ ] 2
- [ ] Part 7 
  - [ ] 1
  - [ ] 2 [DUE TO 05/01]
