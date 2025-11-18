subroutine allocatecpp
    use fields_fkfun
    use conformations
    use chainsdat
    
    implicit none

    ALLOCATE(pxA(cuantasA, nsegA, maxcpp))
    ALLOCATE(pyA(cuantasA, nsegA, maxcpp))
    ALLOCATE(pzA(cuantasA, nsegA, maxcpp))
    
    ALLOCATE(pxB(cuantasB, nsegB, maxcpp))
    ALLOCATE(pyB(cuantasB, nsegB, maxcpp))
    ALLOCATE(pzB(cuantasB, nsegB, maxcpp))

    ALLOCATE(proA(cuantasA, maxcpp))
    ALLOCATE(proB(cuantasB, maxcpp))

    ALLOCATE(zfinalA(cuantasA, maxcpp))
    ALLOCATE(zfinalB(cuantasB, maxcpp))

end subroutine
