!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: simplesampler.F90
!
! !DESCRIPTION: Approximate TROPOMI/Sentinel-5P swath sampling.  Given a UTC
!  hour of day and a (lat,lon) location, returns whether that location falls
!  within the satellite overpass swath.  Used by History_Netcdf_Write to mask
!  out-of-swath grid boxes with NaN before writing to disk.
!\\
!\\
! !INTERFACE:
!
MODULE TropomiSwathModule
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: PointInTropomiSwath
  PUBLIC :: BuildTropomiSwathMask
  PUBLIC :: IsSampledCollection
!
! !PRIVATE TYPES:
!
  REAL(8), PARAMETER :: PI              = 3.14159265358979323846_8
  REAL(8), PARAMETER :: DEG2RAD         = PI / 180.0_8
  REAL(8), PARAMETER :: EARTH_RADIUS_KM = 6371.0_8

  ! --- TROPOMI / Sentinel-5P orbit parameters ---
  REAL(8), PARAMETER :: SWATH_WIDTH_KM        = 2600.0_8
  REAL(8), PARAMETER :: INCLINATION_DEG       = 98.7_8
  REAL(8), PARAMETER :: ASCENDING_NODE_LST_HR = 13.5_8
  REAL(8), PARAMETER :: MARGIN_FACTOR         = 1.6_8

  !-------------------------------------------------------------------------
  ! %%% EDIT THIS LIST %%%
  !
  ! Names of the HISTORY COLLECTIONS to which swath sampling is applied.
  ! Any collection NOT named here is written out in full, unmodified.
  !
  ! DO NOT add 'Restart' to this list -- NaN-ing the restart file would
  ! poison the next simulation segment.
  !-------------------------------------------------------------------------
  INTEGER,           PARAMETER :: N_SAMPLED_COLLECTIONS = 1
  CHARACTER(LEN=63), PARAMETER ::                                            &
       SAMPLED_COLLECTIONS(N_SAMPLED_COLLECTIONS) = (/                       &
          'SpeciesConc                                                    '  &
       /)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: IsSampledCollection
!
! !DESCRIPTION: Returns TRUE if the named HISTORY COLLECTION should have
!  swath sampling applied to it.
!\\
!\\
! !INTERFACE:
!
  FUNCTION IsSampledCollection( Name ) RESULT( IsSampled )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN) :: Name       ! Collection name
!
! !RETURN VALUE:
!
    LOGICAL                      :: IsSampled
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER :: N

    IsSampled = .FALSE.
    DO N = 1, N_SAMPLED_COLLECTIONS
       IF ( TRIM( Name ) == TRIM( SAMPLED_COLLECTIONS(N) ) ) THEN
          IsSampled = .TRUE.
          RETURN
       ENDIF
    ENDDO

  END FUNCTION IsSampledCollection
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: BuildTropomiSwathMask
!
! !DESCRIPTION: Fills a 2-D logical mask (TRUE = in swath) for the given
!  longitude and latitude center arrays at the given UTC hour.  Call this
!  once per file write and reuse the mask for every HISTORY ITEM, rather
!  than calling PointInTropomiSwath box-by-box for each item.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE BuildTropomiSwathMask( HourUTC, Lon, Lat, Mask )
!
! !INPUT PARAMETERS:
!
    REAL(8), INTENT(IN)  :: HourUTC          ! Hour of day, UTC [0,24)
    REAL(8), INTENT(IN)  :: Lon(:)           ! Lon centers [deg east]
    REAL(8), INTENT(IN)  :: Lat(:)           ! Lat centers [deg north]
!
! !OUTPUT PARAMETERS:
!
    LOGICAL, INTENT(OUT) :: Mask(:,:)        ! TRUE = grid box is in swath
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER :: I, J

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J   )
    DO J = 1, SIZE( Lat )
    DO I = 1, SIZE( Lon )
       CALL PointInTropomiSwath( HourUTC, Lat(J), Lon(I), Mask(I,J) )
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE BuildTropomiSwathMask
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: PointInTropomiSwath
!
! !DESCRIPTION: Approximate check for whether (Latitude, Longitude) falls
!  within the TROPOMI/Sentinel-5P swath at the given UTC hour of day.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PointInTropomiSwath( HourUTC, Latitude, Longitude, InSwath )
!
! !INPUT PARAMETERS:
!
    REAL(8), INTENT(IN)  :: HourUTC    ! Hour of day, UTC, in [0,24)
    REAL(8), INTENT(IN)  :: Latitude   ! degrees
    REAL(8), INTENT(IN)  :: Longitude  ! degrees, [-180,180]
!
! !OUTPUT PARAMETERS:
!
    LOGICAL, INTENT(OUT) :: InSwath
!
! !REMARKS:
!  Models the orbit as an idealized great circle inclined relative to the
!  equator, whose ascending-node crossing longitude is derived from the
!  requirement that its local solar time is 13:30 (S5P).  The swath is
!  treated as a band of constant angular half-width straddling that great
!  circle.  Only the ascending-node pass is considered (descending-node
!  pass is intentionally ignored).
!
!  This calculation is intentionally conservative/approximate: no true
!  orbit propagation is performed, and margins are applied so that
!  borderline points are more likely to be included than excluded.
!EOP
!------------------------------------------------------------------------------
!BOC
    REAL(8) :: Lon0Deg, Lon0
    REAL(8) :: Inclination
    REAL(8) :: LatRad, LonRad
    REAL(8) :: Px, Py, Pz
    REAL(8) :: Nx, Ny, Nz
    REAL(8) :: E1x, E1y
    REAL(8) :: DotN
    REAL(8) :: CrossTrackKM
    REAL(8) :: AlongSide
    REAL(8) :: HalfWidthKM

    ! --- Ascending-node longitude implied by the 13:30 LST condition ---
    Lon0Deg = MODULO( 15.0_8 * ( ASCENDING_NODE_LST_HR - HourUTC ), 360.0_8 )
    IF ( Lon0Deg > 180.0_8 ) Lon0Deg = Lon0Deg - 360.0_8
    Lon0 = Lon0Deg * DEG2RAD

    ! Supplementary angle: for a retrograde (inclination > 90 deg)
    ! sun-synchronous orbit, using (180 - inclination) in the tilted-
    ! plane geometry below matches the real ground-track slant direction.
    Inclination = ( 180.0_8 - INCLINATION_DEG ) * DEG2RAD

    ! --- Query point as a unit vector on the sphere ---
    LatRad = Latitude  * DEG2RAD
    LonRad = Longitude * DEG2RAD

    Px = COS(LatRad) * COS(LonRad)
    Py = COS(LatRad) * SIN(LonRad)
    Pz = SIN(LatRad)

    ! --- Orbital-plane pole vector N, tilted by Inclination about the
    !     ascending-node line ---
    Nx = -SIN(Inclination) * SIN(Lon0)
    Ny =  SIN(Inclination) * COS(Lon0)
    Nz =  COS(Inclination)

    ! --- Ascending-node direction vector E1 (in equatorial plane) ---
    E1x = COS(Lon0)
    E1y = SIN(Lon0)

    ! --- Cross-track angular distance: angle between P and the orbital
    !     great circle's plane ---
    DotN = Px * Nx + Py * Ny + Pz * Nz
    DotN = MAX( -1.0_8, MIN( 1.0_8, DotN ) )
    CrossTrackKM = ASIN(DotN) * EARTH_RADIUS_KM

    ! --- Restrict to the ascending arc only (exclude the diametrically
    !     opposite descending arc on the same great circle) ---
    AlongSide = Px * E1x + Py * E1y

    HalfWidthKM = ( SWATH_WIDTH_KM / 2.0_8 ) * MARGIN_FACTOR

    InSwath = ( ABS(CrossTrackKM) <= HalfWidthKM ) .and. ( AlongSide > 0.0_8 )

  END SUBROUTINE PointInTropomiSwath
!EOC
END MODULE TropomiSwathModule
