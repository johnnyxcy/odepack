#ifndef odepack_dlsoda_h
#define odepack_dlsoda_h

#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN
#endif

EXTERN typedef int (*dlsoda_f)(int *neq, double *t, double *y, double *ydot,
                               void *dat);

/**
 * @brief
 *
 * @param f name of subroutine for right-hand side vector f.
 * @param neq number of first order ODEs.
 * @param y array of initial values, of length NEQ.
 * @param t the initial value of the independent variable.
 * @param tout first point where output is desired (.ne. T).
 * @param itol 1 or 2 according as ATOL (below) is a scalar or array.
 * @param rtol relative tolerance parameter (scalar).
 * @param atol absolute tolerance parameter (scalar or array).
                the estimated local error in y(i) will be controlled
                so as to be less than
                    EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
                    EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
                Thus the local error test passes if, in each component,
                either the absolute error is less than ATOL (or ATOL(i)),
                or the relative error is less than RTOL.
                Use RTOL = 0.0 for pure absolute error control, and
                use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
                control.  Caution: actual (global) errors may exceed these
                local tolerances, so choose them conservatively.
 * @param itask 1 for normal computation of output values of y at t = TOUT.
 * @param istate integer flag (input and output).  Set ISTATE = 1.
 * @param iopt 0 to indicate no optional inputs used.
 * @param rwork real work array of length at least:
                22 + NEQ * MAX(16, NEQ + 9).
 * @param lrw declared length of RWORK (in user's dimension).
 * @param iwork integer work array of length at least  20 + NEQ.
 * @param liw declared length of IWORK (in user's dimension).
 * @param jac name of subroutine for Jacobian matrix.
 * @param jt Jacobian type indicator.  Set JT = 2.
 * @return EXTERN

The call sequence parameters used for input only are
    F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, JT,
and those used for both input and output are
    Y, T, ISTATE.
The work arrays RWORK and IWORK are also used for conditional and
optional inputs and optional outputs.  (The term output here refers
to the return from Subroutine DLSODA to the user's calling program.)

The legality of input parameters will be thoroughly checked on the
initial call for the problem, but not checked thereafter unless a
change in input parameters is flagged by ISTATE = 3 on input.

The descriptions of the call arguments are as follows.

F      = the name of the user-supplied subroutine defining the
         ODE system.  The system must be put in the first-order
         form dy/dt = f(t,y), where f is a vector-valued function
         of the scalar t and the vector y.  Subroutine F is to
         compute the function f.  It is to have the form
              SUBROUTINE F (NEQ, T, Y, YDOT)
              DOUBLE PRECISION T, Y(*), YDOT(*)
         where NEQ, T, and Y are input, and the array YDOT = f(t,y)
         is output.  Y and YDOT are arrays of length NEQ.
         Subroutine F should not alter Y(1),...,Y(NEQ).
         F must be declared External in the calling program.

         Subroutine F may access user-defined quantities in
         NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
         (dimensioned in F) and/or Y has length exceeding NEQ(1).
         See the descriptions of NEQ and Y below.

         If quantities computed in the F routine are needed
         externally to DLSODA, an extra call to F should be made
         for this purpose, for consistent and accurate results.
         If only the derivative dy/dt is needed, use DINTDY instead.

NEQ    = the size of the ODE system (number of first order
         ordinary differential equations).  Used only for input.
         NEQ may be decreased, but not increased, during the problem.
         If NEQ is decreased (with ISTATE = 3 on input), the
         remaining components of Y should be left undisturbed, if
         these are to be accessed in F and/or JAC.

         Normally, NEQ is a scalar, and it is generally referred to
         as a scalar in this user interface description.  However,
         NEQ may be an array, with NEQ(1) set to the system size.
         (The DLSODA package accesses only NEQ(1).)  In either case,
         this parameter is passed as the NEQ argument in all calls
         to F and JAC.  Hence, if it is an array, locations
         NEQ(2),... may be used to store other integer data and pass
         it to F and/or JAC.  Subroutines F and/or JAmust include
         NEQ in a Dimension statement in that case.

Y      = a real array for the vector of dependent variables, of
         length NEQ or more.  Used for both input and output on the
         first call (ISTATE = 1), and only for output on other calls.
         On the first call, Y must contain the vector of initial
         values.  On output, Y contains the computed solution vector,
         evaluated at T.  If desired, the Y array may be used
         for other purposes between calls to the solver.

         This array is passed as the Y argument in all calls to
         F and JAC.  Hence its length may exceed NEQ, and locations
         Y(NEQ+1),... may be used to store other real data and
         pass it to F and/or JAC.  (The DLSODA package accesses only
         Y(1),...,Y(NEQ).)

T      = the independent variable.  On input, T is used only on the
         first call, as the initial point of the integration.
         on output, after each call, T is the value at which a
         computed solution Y is evaluated (usually the same as TOUT).
         on an error return, T is the farthest point reached.

TOUT   = the next value of t at which a computed solution is desired.
         Used only for input.

         When starting the problem (ISTATE = 1), TOUT may be equal
         to T for one call, then should .ne. T for the next call.
         For the initial t, an input value of TOUT .ne. T is used
         in order to determine the direction of the integration
         (i.e. the algebraic sign of the step sizes) and the rough
         scale of the problem.  Integration in either direction
         (forward or backward in t) is permitted.

         If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
         the first call (i.e. the first call with TOUT .ne. T).
         Otherwise, TOUT is required on every call.

         If ITASK = 1, 3, or 4, the values of TOUT need not be
         monotone, but a value of TOUT which backs up is limited
         to the current internal T interval, whose endpoints are
         TCUR - HU and TCUR (see optional outputs, below, for
         TCUR and HU).

ITOL   = an indicator for the type of error control.  See
         description below under ATOL.  Used only for input.

RTOL   = a relative error tolerance parameter, either a scalar or
         an array of length NEQ.  See description below under ATOL.
         Input only.

ATOL   = an absolute error tolerance parameter, either a scalar or
         an array of length NEQ.  Input only.

            The input parameters ITOL, RTOL, and ATOL determine
         the error control performed by the solver.  The solver will
         control the vector E = (E(i)) of estimated local errors
         in y, according to an inequality of the form
                     max-norm of ( E(i)/EWT(i) )   .le.   1,
         where EWT = (EWT(i)) is a vector of positive error weights.
         The values of RTOL and ATOL should all be non-negative.
         The following table gives the types (scalar/array) of
         RTOL and ATOL, and the corresponding form of EWT(i).

            ITOL    RTOL       ATOL          EWT(i)
             1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
             2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
             3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
             4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)

         When either of these parameters is a scalar, it need not
         be dimensioned in the user's calling program.

         If none of the above choices (with ITOL, RTOL, and ATOL
         fixed throughout the problem) is suitable, more general
         error controls can be obtained by substituting a
         user-supplied routine for the setting of EWT.
         See Part 4 below.

         If global errors are to be estimated by making a repeated
         run on the same problem with smaller tolerances, then all
         components of RTOL and ATOL (i.e. of EWT) should be scaled
         down uniformly.

ITASK  = an index specifying the task to be performed.
         Input only.  ITASK has the following values and meanings.
         1  means normal computation of output values of y(t) at
            t = TOUT (by overshooting and interpolating).
         2  means take one step only and return.
         3  means stop at the first internal mesh point at or
            beyond t = TOUT and return.
         4  means normal computation of output values of y(t) at
            t = TOUT but without overshooting t = TCRIT.
            TCRIT must be input as RWORK(1).  TCRIT may be equal to
            or beyond TOUT, but not behind it in the direction of
            integration.  This option is useful if the problem
            has a singularity at or beyond t = TCRIT.
         5  means take one step, without passing TCRIT, and return.
            TCRIT must be input as RWORK(1).

         Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
         (within roundoff), it will return T = TCRIT (exactly) to
         indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
         in which case answers at t = TOUT are returned first).

ISTATE = an index used for input and output to specify the
         the state of the calculation.

         On input, the values of ISTATE are as follows.
         1  means this is the first call for the problem
            (initializations will be done).  See note below.
         2  means this is not the first call, and the calculation
            is to continue normally, with no change in any input
            parameters except possibly TOUT and ITASK.
            (If ITOL, RTOL, and/or ATOL are changed between calls
            with ISTATE = 2, the new values will be used but not
            tested for legality.)
         3  means this is not the first call, and the
            calculation is to continue normally, but with
            a change in input parameters other than
            TOUT and ITASK.  Changes are allowed in
            NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, JT, ML, MU,
            and any optional inputs except H0, MXORDN, and MXORDS.
            (See IWORK description for ML and MU.)
         Note:  A preliminary call with TOUT = T is not counted
         as a first call here, as no initialization or checking of
         input is done.  (Such a call is sometimes useful for the
         purpose of outputting the initial conditions.)
         Thus the first call for which TOUT .ne. T requires
         ISTATE = 1 on input.

         On output, ISTATE has the following values and meanings.
          1  means nothing was done; TOUT = T and ISTATE = 1 on input.
          2  means the integration was performed successfully.
         -1  means an excessive amount of work (more than MXSTEP
             steps) was done on this call, before completing the
             requested task, but the integration was otherwise
             successful as far as T.  (MXSTEP is an optional input
             and is normally 500.)  To continue, the user may
             simply reset ISTATE to a value .gt. 1 and call again
             (the excess work step counter will be reset to 0).
             In addition, the user may increase MXSTEP to avoid
             this error return (see below on optional inputs).
         -2  means too much accuracy was requested for the precision
             of the machine being used.  This was detected before
             completing the requested task, but the integration
             was successful as far as T.  To continue, the tolerance
             parameters must be reset, and ISTATE must be set
             to 3.  The optional output TOLSF may be used for this
             purpose.  (Note: If this condition is detected before
             taking any steps, then an illegal input return
             (ISTATE = -3) occurs instead.)
         -3  means illegal input was detected, before taking any
             integration steps.  See written message for details.
             Note:  If the solver detects an infinite loop of calls
             to the solver with illegal input, it will cause
             the run to stop.
         -4  means there were repeated error test failures on
             one attempted step, before completing the requested
             task, but the integration was successful as far as T.
             The problem may have a singularity, or the input
             may be inappropriate.
         -5  means there were repeated convergence test failures on
             one attempted step, before completing the requested
             task, but the integration was successful as far as T.
             This may be caused by an inaccurate Jacobian matrix,
             if one is being used.
         -6  means EWT(i) became zero for some i during the
             integration.  Pure relative error control (ATOL(i)=0.0)
             was requested on a variable which has now vanished.
             The integration was successful as far as T.
         -7  means the length of RWORK and/or IWORK was too small to
             proceed, but the integration was successful as far as T.
             This happens when DLSODA chooses to switch methods
             but LRW and/or LIW is too small for the new method.

         Note:  Since the normal output value of ISTATE is 2,
         it does not need to be reset for normal continuation.
         Also, since a negative input value of ISTATE will be
         regarded as illegal, a negative output value requires the
         user to change it, and possibly other inputs, before
         calling the solver again.

IOPT   = an integer flag to specify whether or not any optional
         inputs are being used on this call.  Input only.
         The optional inputs are listed separately below.
         IOPT = 0 means no optional inputs are being used.
                  default values will be used in all cases.
         IOPT = 1 means one or more optional inputs are being used.

RWORK  = a real array (double precision) for work space, and (in the
         first 20 words) for conditional and optional inputs and
         optional outputs.
         As DLSODA switches automatically between stiff and nonstiff
         methods, the required length of RWORK can change during the
         problem.  Thus the RWORK array passed to DLSODA can either
         have a static (fixed) length large enough for both methods,
         or have a dynamic (changing) length altered by the calling
         program in response to output from DLSODA.

                      --- Fixed Length Case ---
         If the RWORK length is to be fixed, it should be at least
              MAX (LRN, LRS),
         where LRN and LRS are the RWORK lengths required when the
         current method is nonstiff or stiff, respectively.

         The separate RWORK length requirements LRN and LRS are
         as follows:
         IF NEQ is constant and the maximum method orders have
         their default values, then
            LRN = 20 + 16*NEQ,
            LRS = 22 + 9*NEQ + NEQ**2           if JT = 1 or 2,
            LRS = 22 + 10*NEQ + (2*ML+MU)*NEQ   if JT = 4 or 5.
         Under any other conditions, LRN and LRS are given by:
            LRN = 20 + NYH*(MXORDN+1) + 3*NEQ,
            LRS = 20 + NYH*(MXORDS+1) + 3*NEQ + LMAT,
         where
            NYH    = the initial value of NEQ,
            MXORDN = 12, unless a smaller value is given as an
                     optional input,
            MXORDS = 5, unless a smaller value is given as an
                     optional input,
            LMAT   = length of matrix work space:
            LMAT   = NEQ**2 + 2              if JT = 1 or 2,
            LMAT   = (2*ML + MU + 1)*NEQ + 2 if JT = 4 or 5.

                      --- Dynamic Length Case ---
         If the length of RWORK is to be dynamic, then it should
         be at least LRN or LRS, as defined above, depending on the
         current method.  Initially, it must be at least LRN (since
         DLSODA starts with the nonstiff method).  On any return
         from DLSODA, the optional output MCUR indicates the current
         method.  If MCUR differs from the value it had on the
         previous return, or if there has only been one call to
         DLSODA and MCUR is now 2, then DLSODA has switched
         methods during the last call, and the length of RWORK
         should be reset (to LRN if MCUR = 1, or to LRS if
         MCUR = 2).  (An increase in the RWORK length is required
         if DLSODA returned ISTATE = -7, but not otherwise.)
         After resetting the length, call DLSODA with ISTATE = 3
         to signal that change.

LRW    = the length of the array RWORK, as declared by the user.
         (This will be checked by the solver.)

IWORK  = an integer array for work space.
         As DLSODA switches automatically between stiff and nonstiff
         methods, the required length of IWORK can change during
         problem, between
            LIS = 20 + NEQ   and   LIN = 20,
         respectively.  Thus the IWORK array passed to DLSODA can
         either have a fixed length of at least 20 + NEQ, or have a
         dynamic length of at least LIN or LIS, depending on the
         current method.  The comments on dynamic length under
         RWORK above apply here.  Initially, this length need
         only be at least LIN = 20.

         The first few words of IWORK are used for conditional and
         optional inputs and optional outputs.

         The following 2 words in IWORK are conditional inputs:
           IWORK(1) = ML     these are the lower and upper
           IWORK(2) = MU     half-bandwidths, respectively, of the
                      banded Jacobian, excluding the main diagonal.
                      The band is defined by the matrix locations
                      (i,j) with i-ML .le. j .le. i+MU.  ML and MU
                      must satisfy  0 .le.  ML,MU  .le. NEQ-1.
                      These are required if JT is 4 or 5, and
                      ignored otherwise.  ML and MU may in fact be
                      the band parameters for a matrix to which
                      df/dy is only approximately equal.

LIW    = the length of the array IWORK, as declared by the user.
         (This will be checked by the solver.)

Note: The base addresses of the work arrays must not be
altered between calls to DLSODA for the same problem.
The contents of the work arrays must not be altered
between calls, except possibly for the conditional and
optional inputs, and except for the last 3*NEQ words of RWORK.
The latter space is used for internal scratch space, and so is
available for use by the user outside DLSODA between calls, if
desired (but not for use by F or JAC).

JA   = the name of the user-supplied routine to compute the
         Jacobian matrix, df/dy, if JT = 1 or 4.  The JAroutine
         is optional, but if the problem is expected to be stiff much
         of the time, you are encouraged to supply JAC, for the sake
         of efficiency.  (Alternatively, set JT = 2 or 5 to have
         DLSODA compute df/dy internally by difference quotients.)
         If and when DLSODA uses df/dy, it treats this NEQ by NEQ
         matrix either as full (JT = 1 or 2), or as banded (JT =
         4 or 5) with half-bandwidths ML and MU (discussed under
         IWORK above).  In either case, if JT = 1 or 4, the JA
         routine must compute df/dy as a function of the scalar t
         and the vector y.  It is to have the form
              SUBROUTINE JA(NEQ, T, Y, ML, MU, PD, NROWPD)
              DOUBLE PRECISION T, Y(*), PD(NROWPD,*)
         where NEQ, T, Y, ML, MU, and NROWPD are input and the array
         PD is to be loaded with partial derivatives (elements of
         the Jacobian matrix) on output.  PD must be given a first
         dimension of NROWPD.  T and Y have the same meaning as in
         Subroutine F.
              In the full matrix case (JT = 1), ML and MU are
         ignored, and the Jacobian is to be loaded into PD in
         columnwise manner, with df(i)/dy(j) loaded into PD(i,j).
              In the band matrix case (JT = 4), the elements
         within the band are to be loaded into PD in columnwise
         manner, with diagonal lines of df/dy loaded into the rows
         of PD.  Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).
         ML and MU are the half-bandwidth parameters (see IWORK).
         The locations in PD in the two triangular areas which
         correspond to nonexistent matrix elements can be ignored
         or loaded arbitrarily, as they are overwritten by DLSODA.
              JAneed not provide df/dy exactly.  A crude
         approximation (possibly with a smaller bandwidth) will do.
              In either case, PD is preset to zero by the solver,
         so that only the nonzero elements need be loaded by JAC.
         Each call to JAis preceded by a call to F with the same
         arguments NEQ, T, and Y.  Thus to gain some efficiency,
         intermediate quantities shared by both calculations may be
         saved in a user Common block by F and not recomputed by JAC,
         if desired.  Also, JAmay alter the Y array, if desired.
         JAmust be declared External in the calling program.
              Subroutine JAmay access user-defined quantities in
         NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
         (dimensioned in JAC) and/or Y has length exceeding NEQ(1).
         See the descriptions of NEQ and Y above.

JT     = Jacobian type indicator.  Used only for input.
         JT specifies how the Jacobian matrix df/dy will be
         treated, if and when DLSODA requires this matrix.
         JT has the following values and meanings:
          1 means a user-supplied full (NEQ by NEQ) Jacobian.
          2 means an internally generated (difference quotient) full
            Jacobian (using NEQ extra calls to F per df/dy value).
          4 means a user-supplied banded Jacobian.
          5 means an internally generated banded Jacobian (using
            ML+MU+1 extra calls to F per df/dy evaluation).
         If JT = 1 or 4, the user must supply a Subroutine JA
         (the name is arbitrary) as described above under JAC.
         If JT = 2 or 5, a dummy argument can be used.
-----------------------------------------------------------------------
Optional Inputs.

The following is a list of the optional inputs provided for in the
call sequence.  (See also Part 2.)  For each such input variable,
this table lists its name as used in this documentation, its
location in the call sequence, its meaning, and the default value.
The use of any of these inputs requires IOPT = 1, and in that
case all of these inputs are examined.  A value of zero for any
of these optional inputs will cause the default value to be used.
Thus to use a subset of the optional inputs, simply preload
locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
then set those of interest to nonzero values.

Name    Location      Meaning and Default Value

H0      RWORK(5)  the step size to be attempted on the first step.
                  The default value is determined by the solver.

HMAX    RWORK(6)  the maximum absolute step size allowed.
                  The default value is infinite.

HMIN    RWORK(7)  the minimum absolute step size allowed.
                  The default value is 0.  (This lower bound is not
                  enforced on the final step before reaching TCRIT
                  when ITASK = 4 or 5.)

IXPR    IWORK(5)  flag to generate extra printing at method switches.
                  IXPR = 0 means no extra printing (the default).
                  IXPR = 1 means print data on each switch.
                  T, H, and NST will be printed on the same logical
                  unit as used for error messages.

MXSTEP  IWORK(6)  maximum number of (internally defined) steps
                  allowed during one call to the solver.
                  The default value is 500.

MXHNIL  IWORK(7)  maximum number of messages printed (per problem)
                  warning that T + H = T on a step (H = step size).
                  This must be positive to result in a non-default
                  value.  The default value is 10.

MXORDN  IWORK(8)  the maximum order to be allowed for the nonstiff
                  (Adams) method.  the default value is 12.
                  if MXORDN exceeds the default value, it will
                  be reduced to the default value.
                  MXORDN is held constant during the problem.

MXORDS  IWORK(9)  the maximum order to be allowed for the stiff
                  (BDF) method.  The default value is 5.
                  If MXORDS exceeds the default value, it will
                  be reduced to the default value.
                  MXORDS is held constant during the problem.
-----------------------------------------------------------------------
Optional Outputs.

As optional additional output from DLSODA, the variables listed
below are quantities related to the performance of DLSODA
which are available to the user.  These are communicated by way of
the work arrays, but also have internal mnemonic names as shown.
except where stated otherwise, all of these outputs are defined
on any successful return from DLSODA, and on any return with
ISTATE = -1, -2, -4, -5, or -6.  On an illegal input return
(ISTATE = -3), they will be unchanged from their existing values
(if any), except possibly for TOLSF, LENRW, and LENIW.
On any error return, outputs relevant to the error will be defined,
as noted below.

Name    Location      Meaning

HU      RWORK(11) the step size in t last used (successfully).

HCUR    RWORK(12) the step size to be attempted on the next step.

TCUR    RWORK(13) the current value of the independent variable
                  which the solver has actually reached, i.e. the
                  current internal mesh point in t.  On output, TCUR
                  will always be at least as far as the argument
                  T, but may be farther (if interpolation was done).

TOLSF   RWORK(14) a tolerance scale factor, greater than 1.0,
                  computed when a request for too much accuracy was
                  detected (ISTATE = -3 if detected at the start of
                  the problem, ISTATE = -2 otherwise).  If ITOL is
                  left unaltered but RTOL and ATOL are uniformly
                  scaled up by a factor of TOLSF for the next call,
                  then the solver is deemed likely to succeed.
                  (The user may also ignore TOLSF and alter the
                  tolerance parameters in any other way appropriate.)

TSW     RWORK(15) the value of t at the time of the last method
                  switch, if any.

NST     IWORK(11) the number of steps taken for the problem so far.

NFE     IWORK(12) the number of f evaluations for the problem so far.

NJE     IWORK(13) the number of Jacobian evaluations (and of matrix
                  LU decompositions) for the problem so far.

NQU     IWORK(14) the method order last used (successfully).

NQCUR   IWORK(15) the order to be attempted on the next step.

IMXER   IWORK(16) the index of the component of largest magnitude in
                  the weighted local error vector ( E(i)/EWT(i) ),
                  on an error return with ISTATE = -4 or -5.

LENRW   IWORK(17) the length of RWORK actually required, assuming
                  that the length of RWORK is to be fixed for the
                  rest of the problem, and that switching may occur.
                  This is defined on normal returns and on an illegal
                  input return for insufficient storage.

LENIW   IWORK(18) the length of IWORK actually required, assuming
                  that the length of IWORK is to be fixed for the
                  rest of the problem, and that switching may occur.
                  This is defined on normal returns and on an illegal
                  input return for insufficient storage.

MUSED   IWORK(19) the method indicator for the last successful step:
                  1 means Adams (nonstiff), 2 means BDF (stiff).

MCUR    IWORK(20) the current method indicator:
                  1 means Adams (nonstiff), 2 means BDF (stiff).
                  This is the method to be attempted
                  on the next step.  Thus it differs from MUSED
                  only if a method switch has just been made.

The following two arrays are segments of the RWORK array which
may also be of interest to the user as optional outputs.
For each array, the table below gives its internal name,
its base address in RWORK, and its description.

Name    Base Address      Description

YH      21             the Nordsieck history array, of size NYH by
                       (NQCUR + 1), where NYH is the initial value
                       of NEQ.  For j = 0,1,...,NQCUR, column j+1
                       of YH contains HCUR**j/factorial(j) times
                       the j-th derivative of the interpolating
                       polynomial currently representing the solution,
                       evaluated at T = TCUR.

ACOR     LACOR         array of size NEQ used for the accumulated
        (from Common   corrections on each step, scaled on output
          as noted)    to represent the estimated local error in y
                       on the last step.  This is the vector E in
                       the description of the error control.  It is
                       defined only on a successful return from
                       DLSODA.  The base address LACOR is obtained by
                       including in the user's program the
                       following 2 lines:
                          COMMON /DLS001/ RLS(218), ILS(37)
                          LACOR = ILS(22)

-----------------------------------------------------------------------
Part 2.  Other Routines Callable.

The following are optional calls which the user may make to
gain additional capabilities in conjunction with DLSODA.
(The routines XSETUN and XSETF are designed to conform to the
SLATEerror handling package.)

    Form of Call                  Function
  CALL XSETUN(LUN)          set the logical unit number, LUN, for
                            output of messages from DLSODA, if
                            the default is not desired.
                            The default value of LUN is 6.

  CALL XSETF(MFLAG)         set a flag to control the printing of
                            messages by DLSODA.
                            MFLAG = 0 means do not print. (Danger:
                            This risks losing valuable information.)
                            MFLAG = 1 means print (the default).

                            Either of the above calls may be made at
                            any time and will take effect immediately.

  CALL DSRCMA(RSAV,ISAV,JOB) saves and restores the contents of
                            the internal Common blocks used by
                            DLSODA (see Part 3 below).
                            RSAV must be a real array of length 240
                            or more, and ISAV must be an integer
                            array of length 46 or more.
                            JOB=1 means save Common into RSAV/ISAV.
                            JOB=2 means restore Common from RSAV/ISAV.
                               DSRCMA is useful if one is
                            interrupting a run and restarting
                            later, or alternating between two or
                            more problems solved with DLSODA.

  CALL DINTDY(,,,,,)        provide derivatives of y, of various
       (see below)          orders, at a specified point t, if
                            desired.  It may be called only after
                            a successful return from DLSODA.

The detailed instructions for using DINTDY are as follows.
The form of the call is:

  CALL DINTDY (T, K, RWORK(21), NYH, DKY, IFLAG)

The input parameters are:

T         = value of independent variable where answers are desired
            (normally the same as the T last returned by DLSODA).
            For valid results, T must lie between TCUR - HU and TCUR.
            (See optional outputs for TCUR and HU.)
K         = integer order of the derivative desired.  K must satisfy
            0 .le. K .le. NQCUR, where NQCUR is the current order
            (see optional outputs).  The capability corresponding
            to K = 0, i.e. computing y(T), is already provided
            by DLSODA directly.  Since NQCUR .ge. 1, the first
            derivative dy/dt is always available with DINTDY.
RWORK(21) = the base address of the history array YH.
NYH       = column length of YH, equal to the initial value of NEQ.

The output parameters are:

DKY       = a real array of length NEQ containing the computed value
            of the K-th derivative of y(t).
IFLAG     = integer flag, returned as 0 if K and T were legal,
            -1 if K was illegal, and -2 if T was illegal.
            On an error return, a message is also written.
-----------------------------------------------------------------------
Part 3.  Common Blocks.

If DLSODA is to be used in an overlay situation, the user
must declare, in the primary overlay, the variables in:
  (1) the call sequence to DLSODA, and
  (2) the two internal Common blocks
        /DLS001/  of length  255  (218 double precision words
                     followed by 37 integer words),
        /DLSA01/  of length  31    (22 double precision words
                     followed by  9 integer words).

If DLSODA is used on a system in which the contents of internal
Common blocks are not preserved between calls, the user should
declare the above Common blocks in the calling program to insure
that their contents are preserved.

If the solution of a given problem by DLSODA is to be interrupted
and then later continued, such as when restarting an interrupted run
or alternating between two or more problems, the user should save,
following the return from the last DLSODA call prior to the
interruption, the contents of the call sequence variables and the
internal Common blocks, and later restore these values before the
next DLSODA call for that problem.  To save and restore the Common
blocks, use Subroutine DSRCMA (see Part 2 above).

-----------------------------------------------------------------------
Part 4.  Optionally Replaceable Solver Routines.

Below is a description of a routine in the DLSODA package which
relates to the measurement of errors, and can be
replaced by a user-supplied version, if desired.  However, since such
a replacement may have a major impact on performance, it should be
done only when absolutely necessary, and only with great caution.
(Note: The means by which the package version of a routine is
superseded by the user's version may be system-dependent.)

(a) DEWSET.
The following subroutine is called just before each internal
integration step, and sets the array of error weights, EWT, as
described under ITOL/RTOL/ATOL above:
    Subroutine DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
where NEQ, ITOL, RTOL, and ATOL are as in the DLSODA call sequence,
YCUR contains the current dependent variable vector, and
EWT is the array of weights set by DEWSET.

If the user supplies this subroutine, it must return in EWT(i)
(i = 1,...,NEQ) a positive quantity suitable for comparing errors
in y(i) to.  The EWT array returned by DEWSET is passed to the
DMNORM routine, and also used by DLSODA in the computation
of the optional output IMXER, and the increments for difference
quotient Jacobians.

In the user-supplied version of DEWSET, it may be desirable to use
the current values of derivatives of y.  Derivatives up to order NQ
are available from the history array YH, described above under
optional outputs.  In DEWSET, YH is identical to the YCUR array,
extended to NQ + 1 columns with a column length of NYH and scale
factors of H**j/factorial(j).  On the first call for the problem,
given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
NYH is the initial value of NEQ.  The quantities NQ, H, and NST
can be obtained by including in DEWSET the statements:
    DOUBLE PRECISION RLS
    COMMON /DLS001/ RLS(218),ILS(37)
    NQ = ILS(33)
    NST = ILS(34)
    H = RLS(212)
Thus, for example, the current value of dy/dt can be obtained as
YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
unnecessary when NST = 0).
-----------------------------------------------------------------------

C***REVISION HISTORY  (YYYYMMDD)
19811102  DATE WRITTEN
19820126  Fixed bug in tests of work space lengths;
          minor corrections in main prologue and comments.
19870330  Major update: corrected comments throughout;
          removed TRET from Common; rewrote EWSET with 4 loops;
          fixed t test in INTDY; added Cray directives in STODA;
          in STODA, fixed DELP init. and logic around PJAcall;
          combined routines to save/restore Common;
          passed LEVEL = 0 in error message calls (except run abort).
19970225  Fixed lines setting JSTART = -2 in Subroutine LSODA.
20010425  Major update: convert source lines to upper case;
          added *DECK lines; changed from 1 to * in dummy dimensions;
          changed names R1MACH/D1MACH to RUMACH/DUMACH;
          renamed routines for uniqueness across single/double prec.;
          converted intrinsic names to generic form;
          removed ILLIN and NTREP (data loaded) from Common;
          removed all 'own' variables from Common;
          changed error messages to quoted strings;
          replaced XERRWV/XERRWD with 1993 revised version;
          converted prologues, comments, error messages to mixed case;
          numerous corrections to prologues and internal comments.
20010507  Converted single precision source to double precision.
20010613  Revised excess accuracy test (to match rest of ODEPACK).
20010808  Fixed bug in DPRJA (matrix in DBNORM call).
20020502  Corrected declarations in descriptions of user routines.
20031105  Restored 'own' variables to Common blocks, to enable
          interrupt/restart feature.
20031112  Added SAVE statements for data-loaded constants.

-----------------------------------------------------------------------
Other routines in the DLSODA package.

In addition to Subroutine DLSODA, the DLSODA package includes the
following subroutines and function routines:
 DINTDY   computes an interpolated value of the y vector at t = TOUT.
 DSTODA   is the core integrator, which does one step of the
          integration and the associated error control.
 DCFODE   sets all method coefficients and test constants.
 DPRJA    computes and preprocesses the Jacobian matrix J = df/dy
          and the Newton iteration matrix P = I - h*l0*J.
 DSOLSY   manages solution of linear system in chord iteration.
 DEWSET   sets the error weight vector EWT before each step.
 DMNORM   computes the weighted max-norm of a vector.
 DFNORM   computes the norm of a full matrix consistent with the
          weighted max-norm on vectors.
 DBNORM   computes the norm of a band matrix consistent with the
          weighted max-norm on vectors.
 DSRCMA   is a user-callable routine to save and restore
          the contents of the internal Common blocks.
 DGEFA and DGESL   are routines from LINPACK for solving full
          systems of linear algebraic equations.
 DGBFA and DGBSL   are routines from LINPACK for solving banded
          linear systems.
 DUMACH   computes the unit roundoff in a machine-independent manner.
 XERRWD, XSETUN, XSETF, IXSAV, and IUMACH  handle the printing of all
          error messages and warnings.  XERRWD is machine-dependent.
Note:  DMNORM, DFNORM, DBNORM, DUMACH, IXSAV, and IUMACH are
function routines.  All the others are subroutines.

-----------------------------------------------------------------------
 */
EXTERN void dlsoda(dlsoda_f f, int *neq, double *y, double *t, double *tout,
                   int *itol, double *rtol, double *atol, int *itask,
                   int *istate, int *iopt, double *rwork, int *lrw, int *iwork,
                   int *liw,
                   void (*jac)(int *neq, double *t, double *y, int *ml, int *mu,
                               double *pd, int *nrowpd),
                   int *jt, void *dat);

#endif /* odepack_dlsoda_h */