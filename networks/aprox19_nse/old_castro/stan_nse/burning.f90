module burning_module

  use bl_types
  use bl_constants_module
  use bl_error_module


contains

   subroutine burningtable_init() 
   use bl_error_module
   use parallel
   implicit none



   double precision rhozerol(73), t9(73), rhof(73), helium(73), &
          carbon(73), oxygen(73), xneon(73), xmagnes(73), silicon(73), & 
          xiron(73), drhodx12(73), bea(73), abar(73), atwood(73)

   common /table_burning/ rhozerol, t9, rhof, helium, &
          carbon, oxygen, xneon, xmagnes, silicon, &
          xiron, drhodx12, bea, abar, atwood 
   save /table_burning/

   integer          :: nden,i

   logical          first1
   data             first1/.false./

   if ( first1 ) return

   first1 = .true.

      
!     set table parameters

      nden=73
    
!     read in table 
      open(unit=4,file='burning_table.dat',status='old',err=1010)
      GO TO 1011
 1010 continue
      call bl_error('EOSINIT: Failed to open burning_table.dat') 
 1011 continue 

      do i = 1,nden
      read (4,10) rhozerol(i), t9(i), rhof(i), helium(i), carbon(i),& 
          oxygen(i), xneon(i), xmagnes(i), silicon(i), xiron(i),&
          drhodx12(i), bea(i), abar(i), atwood(i)
10 format (1pe12.3,7e12.3)
      enddo

      close(4)

   if ( parallel_ioprocessor() ) then
        write(6,*)
        write(6,*) 'finished reading burning table'
        write(6,*)
   end if

   end subroutine burningtable_init

   subroutine networkburn(rholog,c12,rhologf,hef,c12f,o16f,xnef,xmgf,xsif,xfef,beaf,abarf,atwoodf)

   use bl_error_module
   use parallel
   implicit none

! if this is made into a subroutine, the argments are
!    rholog = log base 10 of current density
!    c12    = current carbon mass fraction
!    hef    = returned final helium mass fraction
!    c12f   = returned final carbon mass fraction in ash
!    o16f   = returned final oxygen mass fraction in ash
!    xnef   = returned final neon mass fraction in ash
!    xmgf   = returned final magnesium mass fraction in ash
!    xsif   = returned final silicon - calcium mass fraction in ash
!    xfef   = returned final iron mass fraction in ash
!    beaf   = returned final binding energy in ash
!    abarf  = returned final abar in ash

   double precision rholog,c12,rhologf,hef,c12f,o16f,xnef,xmgf,xsif,xfef,beaf,abarf,atwoodf

   double precision rhozerol(73), t9(73), rhof(73), helium(73), &
          carbon(73), oxygen(73), xneon(73), xmagnes(73), silicon(73), & 
          xiron(73), drhodx12(73), bea(73), abar(73), atwood(73)

   common /table_burning/ rhozerol, t9, rhof, helium, &
          carbon, oxygen, xneon, xmagnes, silicon, &
          xiron, drhodx12, bea, abar, atwood
   save /table_burning/

   double precision, parameter :: delrho = 0.05d0 
   double precision, parameter :: rhomin = 6.40d0
   double precision            :: drho, drhodx, rhoiter
   integer                     :: i,j

   if (c12.gt.0.50d0) c12 = 0.50d0
   if (c12.lt.0.0d0)  c12 = 0.0d0

   if (rholog.gt.10.0d0) rholog = 10.0d0

   if (rholog.lt.rhomin) then
        hef  = 0.0d0
        c12f = 0.50d0
        o16f = 0.50d0
        xnef = 0.0d0
        xmgf = 0.0d0
        xsif = 0.0d0
        xfef = 0.0d0
        beaf = 7.8283d0
        abarf= 13.7d0
        atwoodf = 0.d0

   else

!   determine what initial density to use based on current c12 and rho
!      first find value for drho/dx12 at this density. Iterate 5 times.
!     More iteration improves accuracy but slows code

      ! Initialize the value of rhoiter before we start the i loop
      rhoiter = rholog

      do  i = 1,5
         j = (rhoiter - rhomin)/delrho + 1.d-10
         j = j+1
         drho = (rhoiter - rhozerol(j))/delrho
         drhodx = drhodx12(j)+ drho*(drhodx12(j+1)-drhodx12(j))
         rhoiter = 10.0d0**(rholog) + drhodx*(0.50d0 - c12)
         rhoiter = log10(rhoiter)
      enddo
      rhologf = rhoiter

!  determine characteristic of ash for this fuel density

   j = (rhologf - rhomin)/delrho + 1.d-10
   j = j+1   
   drho = (rhologf - rhozerol(j))/delrho
   hef = helium(j) + drho *  (helium(j+1)-helium(j))
   c12f = carbon(j) + drho * (carbon(j+1)-carbon(j))
   o16f = oxygen(j) + drho * (oxygen(j+1)-oxygen(j))
   xnef = xneon(j) + drho * (xneon(j+1)-xneon(j))
   xmgf = xmagnes(j) + drho * (xmagnes(j+1)-xmagnes(j))
   xsif = silicon(j) + drho * (silicon(j+1)-silicon(j))
   xfef = xiron(j) + drho * (xiron(j+1)-xiron(j))
   beaf = bea(j) + drho * (bea(j+1)-bea(j))
   abarf = abar(j) + drho * (abar(j+1)-abar(j))
   atwoodf=atwood(j)+drho * (atwood(j+1)-atwood(j)) 

   endif
   end subroutine networkburn

end module burning_module
