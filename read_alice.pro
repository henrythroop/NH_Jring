pro read_alice

; # I am trying to use IDL's FITS reader to read the 'analog' Alice data from HK extension, which Randy 
; insists is there. I can't find it, using either python or IDL. But here is a shot...

common units

file =  '/Users/throop/Data/NH_Alice_Ring/O_RING_OC3/data/pluto/level2/ali/all/ali_0299391413_0x4b5_sci_1.fit'

a5 = mrdfits(file, 5)  

print a5

stop



