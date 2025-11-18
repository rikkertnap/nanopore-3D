subroutine allocatencha

    use fields_fkfun
    use chainsdat
    use conformations
    use rotchain
    implicit none

    print*,"allocatencha:",ncha

    ! fields_fkfun
    !ALLOCATE(q(ncha))
    ALLOCATE(qA(ncha))
    ALLOCATE(sumgaucheA(ncha))
    ALLOCATE(ngaucheA(cuantasA,ncha))
    ALLOCATE(qB(ncha))
    ALLOCATE(sumgaucheB(ncha))
    ALLOCATE(ngaucheB(cuantasB,ncha))

    ! chainsdat
    allocate(posicionA(ncha,3))
    allocate(posicionB(ncha,3))
    
    allocate(ngpol(ncha))

    allocate(newcuantasA(ncha))
    allocate(newcuantasB(ncha))

    allocate(hasGraftA(ncha))    ! = true if graftpoint i has a A type polymer
    allocate(hasGraftB(ncha)) 

end subroutine
