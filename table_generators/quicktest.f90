program quicktest
    use polytools
    use astro_constants

    implicit none
    double precision :: ns,Ms,Bns,beta_c
    ns = 3.D0
    Ms = 2.D33
    Bns = 0.17D0


    call beta_func_solver(ns,Ms,Bns,1.D0,beta_c)


    write(*,*) beta_c

end program quicktest
