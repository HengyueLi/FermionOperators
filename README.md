# FermionOperators
This module is used for building a package for numerical calculation on the spin-1/2 fermion system. For a fermion system, a general state |X> can be expanded by the [basis](https://github.com/HengyueLi/Fermion_Table#representation-of-manybody-state) {|ѱ<sub>i</sub>>}:</br>
|X> = Σ<sub>i</sub> |ѱ<sub>i</sub>>(1).</br>&nbsp;
A general fermion operator f can act on this state and a new state |Y> is obtained:</br>
|Y> = f |X>    (2)  </br>&nbsp;
 Here we define the concept of 'single operator' is a multiplication of c and c<sup>+</sup>, where c and c<sup>+</sup> are the creation and annihilation operators. For example c<sup>+</sup>,c,c<sup>+</sup>c,cccc... are all single operators and (c + c<sup>+</sup>) is not. We can futher expand f as a summation of 'single operators'.</br>
 f = Σ<sub>i</sub>A<sub>i</sub>.</br>
 Here A<sub>i</sub> the single operator and is what this module simulated. We can expand (2) by expression similar to (1):</br>
|Y> = Σ<sub>i</sub> |ѱ'<sub>i</sub>> = Σ<sub>a</sub>Σ<sub>i</sub> A<sub>a</sub>  |ѱ<sub>i</sub>>.</br>
We can conclude all we need is to calculate A<sub>a</sub>  |ѱ<sub>i</sub>>. </br>
In this module, all the usefull actions like: c, c<sup>+</sup>,c<sup>+</sup>c,c<sup>+</sup>cc<sup>+</sup>c...will be offered.
## example:

    program main
      use FermionOperators

      type(FermOper):: f
      integer       :: ns,optid,para(8)
      integer*8     :: psi1,psi2
      integer       :: sign


      ! In this example, we set the size of system to be 3.
      ns    = 3


      ! In this example, we test operator   c^+_i c_j
      ! The details of how to chose optid can be found in source code directly.
      optid = 3




      ! notice that the index of site is marked from 0 to 2 for this system ( ns = 3)
      ! In c^+_i c_j, we set i = { siteid = 0   spin = up  }
      para(1) = 0  ;  para(2) = 0
      ! In c^+_i c_j, we set j = { siteid = 1   spin = down}
      para(3) = 1  ;  para(4) = 1


      ! Initialization :
      call f%Initialization(ns,optid,para)


      ! We set psi1 = 50 for test
      ! |binary(50)> = |110010> = |010011(←)> = |0↑0>⊗|0↓↓>
      ! (see the meanning of '←' in https://github.com/HengyueLi/Fermion_Table)
      psi1 = 50_8


      ! we would expect that after f is act on psi1, the output state should be:
      ! psi2 = -|↑↑0>⊗|00↓> = -|110001(←)> = -|100011>=-|Octal(35)>
      ! notice that there will be a sign=-1 appear

      !  try it
      call f%act(psi1,psi2,sign)



      ! the output should be:
      ! psi2 = 35 , sign = -1   as we expected.
      write(*,*)"result is:"
      write(*,*)"psi2=",psi2
      write(*,*)"the possible sign is=",sign


    endprogram



The following table will be frequently used.
    !  ╔═════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╗
    !  ║  optid = -1:                                                                                                                    ║
    !  ║              always returen unavalabel state   ( output = -1)                                                                   ║
    !  ║  optid = 0 :                                                                                                                    ║
    !  ║              do nothing                                                                                                         ║
    !  ║  optid = 1 : c_{i,spin}                                                                                                         ║
    !  ║              i=para(1) ; spin=para(2)                                                                                           ║
    !  ║                                                                                                                                 ║
    !  ║  optid = 2 : c^+_{i,spin}                                                                                                       ║
    !  ║              i=para(1) ; spin=para(2)                                                                                           ║
    !  ║                                                                                                                                 ║
    !  ║  optid = 3 : c^+_{i,sini} * c_{j,spinj}                                                                                         ║
    !  ║              i=para(1) ; spini=para(2) ; j=para(3) ; spinj=para(4)                                                              ║
    !  ║                                                                                                                                 ║
    !  ║  optid = 4 : c_{i,sini} * c_{j,spinj}                                                                                           ║
    !  ║              i=para(1) ; spini=para(2) ; j=para(3) ; spinj=para(4)                                                              ║
    !  ║                                                                                                                                 ║
    !  ║  optid = 5 : c^+_{i,sini} * c^+_{j,spinj}                                                                                       ║
    !  ║              i=para(1) ; spini=para(2) ; j=para(3) ; spinj=para(4)                                                              ║
    !  ║                                                                                                                                 ║
    !  ║  optid = 6 : c^+_{i,sini} * c_{j,spinj} * c^+_{k,sini} * c_{l,spinj}                                                            ║
    !  ║              i=para(1) ; spini=para(2) ; j=para(3) ; spinj=para(4) ; k=para(5) ; spink=para(6) ; l=para(7) ; spinl=para(8) ;    ║
    !  ║                                                                                                                                 ║
    !  ║  optid = 7 : c^+_{i,sini} * c^+_{j,spinj} * c_{k,sini} * c_{l,spinj}                                                            ║
    !  ║              i=para(1) ; spini=para(2) ; j=para(3) ; spinj=para(4) ; k=para(5) ; spink=para(6) ; l=para(7) ; spinl=para(8) ;    ║
    !  ║                                                                                                                                 ║
    !  ║  optid = 8 : n_{i,spini} * n_{j,spinj}                                                                                          ║
    !  ║              i=para(1) ; spini=para(2) ; j=para(3) ; spinj=para(4)                                                              ║
    !  ║                                                                                                                                 ║
    !  ║  optid = 9 : n_{i,spini}                                                                                                        ║
    !  ║              i=para(1) ; spini=para(2)                                                                                          ║
    !  ╚═════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╝
    !  
