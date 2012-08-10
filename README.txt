RKLab provides a suite of Matlab files that implement adaptive-step
explicit and implicit Runge-Kutta methods.  This software is designed
for instructional use and has been implemented to be easily understood
and modified.  As a result, certain efficiency-related optimizations
have been omitted.

We currently provide three classes of methods:
(a) explicit Runge-Kutta (ERK), 
(b) diagonally-implicit Runge-Kutta (DIRK), and
(c) fully implicit Runge-Kutta (IRK).

In addition to these three classes of methods, we provide a large
number of Butcher tables (over 70), holding the coefficients for
existing methods in each of these three categories.

The two most-developed solver classes are the ERK and DIRK methods.
Most of the included Butcher tables for these methods include
embeddings, allowing for straight-forward time step adaptivity
controls.  However, all three classes of methods may be run in
so-called "fixed-step mode", wherein the solver will take steps
of a user-supplied magnitude.

In addition to these sets of solvers, we provide three example
problems that may be used to test different methods, and that may be
used as a template for creating new problems.  These test problems
include:
(a) a non-stiff variant of the Van der Pol oscillator (two-component
    nonlinear ODE system). 
(b) a stiff Brusselator problem (three-component nonlinear ODE
    system), and
(c) a PDE variant of the Brusselator problem (two-component nonlinear
    PDE system).

