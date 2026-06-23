
C *** INPUT ECOSYSTEM FILES ***

c 1KM ECOSYSTEM 
c --- Fortran OPEN 
      integer		ECO_1KM_UNIT
      parameter		(ECO_1KM_UNIT = 10) 
      character*3	ECO_1KM_STATUS
      parameter		(ECO_1KM_STATUS = 'OLD')
      character*11	ECO_1KM_FORM
      parameter		(ECO_1KM_FORM = 'UNFORMATTED')
      character*6	ECO_1KM_ACCESS
      parameter		(ECO_1KM_ACCESS = 'DIRECT')
      integer		ECO_1KM_RECL
      parameter		(ECO_1KM_RECL = 512)

C *** INPUT DATA FILES ***

c ANC THRESHOLDS
c --- Fortran OPEN
      integer		ANC_THRESHOLD_UNIT     
      parameter		(ANC_THRESHOLD_UNIT = 41)
      character*3	ANC_THRESHOLD_STATUS     
      parameter		(ANC_THRESHOLD_STATUS = 'OLD')
      character*9	ANC_THRESHOLD_FORM     
      parameter		(ANC_THRESHOLD_FORM = 'FORMATTED')
      character*10	ANC_THRESHOLD_ACCESS     
      parameter		(ANC_THRESHOLD_ACCESS = 'SEQUENTIAL')

c ANC OCEAN 0.86 CLOUD TEST THRESHOLDS
c --- Fortran OPEN
      integer		ANC_OCNTHR_UNIT     
      parameter		(ANC_OCNTHR_UNIT = 42)
      character*3	ANC_OCNTHR_STATUS     
      parameter		(ANC_OCNTHR_STATUS = 'OLD')
      character*9	ANC_OCNTHR_FORM     
      parameter		(ANC_OCNTHR_FORM = 'FORMATTED')
      character*10	ANC_OCNTHR_ACCESS     
      parameter		(ANC_OCNTHR_ACCESS = 'SEQUENTIAL')

C *** OUTPUT FILES ***

c RUNTIME DEBUG 
c --- Fortran OPEN 
      integer		DEBUG_UNIT
      parameter		(DEBUG_UNIT = 90) 
      character*7	DEBUG_STATUS
      parameter		(DEBUG_STATUS = 'UNKNOWN')
      character*9	DEBUG_FORM
      parameter		(DEBUG_FORM = 'FORMATTED')
      character*10	DEBUG_ACCESS
      parameter		(DEBUG_ACCESS = 'SEQUENTIAL')

c MOD35 CLOUDMASK  
c --- Fortran OPEN 
      integer		MOD35_MASK_UNIT
      parameter		(MOD35_MASK_UNIT = 91) 
      character*7	MOD35_MASK_STATUS
      parameter		(MOD35_MASK_STATUS = 'UNKNOWN')
      character*11	MOD35_MASK_FORM
      parameter		(MOD35_MASK_FORM = 'UNFORMATTED')
      character*6	MOD35_MASK_ACCESS
      parameter		(MOD35_MASK_ACCESS = 'DIRECT')

c --- Fortran OPEN 
      integer		MOD35_HDR_UNIT
      parameter		(MOD35_HDR_UNIT = 92) 
      character*7	MOD35_HDR_STATUS
      parameter		(MOD35_HDR_STATUS = 'UNKNOWN')
      character*9	MOD35_HDR_FORM
      parameter		(MOD35_HDR_FORM = 'FORMATTED')
      character*10	MOD35_HDR_ACCESS
      parameter		(MOD35_HDR_ACCESS = 'SEQUENTIAL')

c MOD35 QUALITY CONTROL  
c --- Fortran OPEN 
      integer		MOD35_QA_UNIT
      parameter		(MOD35_QA_UNIT = 93) 
      character*7	MOD35_QA_STATUS
      parameter		(MOD35_QA_STATUS = 'UNKNOWN')
      character*11	MOD35_QA_FORM
      parameter		(MOD35_QA_FORM = 'UNFORMATTED')
      character*6	MOD35_QA_ACCESS
      parameter		(MOD35_QA_ACCESS = 'DIRECT')

c --- Fortran OPEN 
      integer		MOD35_QAHDR_UNIT
      parameter		(MOD35_QAHDR_UNIT = 94) 
      character*7	MOD35_QAHDR_STATUS
      parameter		(MOD35_QAHDR_STATUS = 'UNKNOWN')
      character*9	MOD35_QAHDR_FORM
      parameter		(MOD35_QAHDR_FORM = 'FORMATTED')
      character*10	MOD35_QAHDR_ACCESS
      parameter		(MOD35_QAHDR_ACCESS = 'SEQUENTIAL')

