# Fit-Student-t-distribution
Fit a Student t distribution to univariate data via maximum likelihood. The Nelder-Mead algorithm is used for optimization. The degrees-of-freedom parameter of the Student t distribution can either be estimated or fixed. Compile with

`gfortran kind.f90 constants.f90 random.f90 basic_stats.f90 obj_fun_student_t.f90 nelder_mead.f90 student_t_fit.f90 xstudent_t_fit.f90`

```
Sample output:

 #obs:      100000
           parameter:           mu            s           nu         mean
                true:     0.500000     1.200000     3.000000
           estimated:     0.477349     1.212073     3.190233     0.490521
           estimated:     0.497564     0.860047     1.000000     0.490521
           estimated:     0.497169     1.081049     2.000000     0.490521
           estimated:     0.496631     1.202713     3.000000     0.490521
           estimated:     0.499521     1.286711     4.000000     0.490521
           estimated:     0.498273     1.351200     5.000000     0.490521

           estimated:     0.487178     1.200297     2.975005     0.498211
           estimated:     0.491936     0.864333     1.000000     0.498211
           estimated:     0.492263     1.072030     2.000000     0.498211
           estimated:     0.493383     1.201701     3.000000     0.498211
           estimated:     0.500746     1.279792     4.000000     0.498211
           estimated:     0.495067     1.350064     5.000000     0.498211
```
