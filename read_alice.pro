pro read_alice

; # I am trying to use IDL's FITS reader to read the 'analog' Alice data from HK extension, which Randy 
; insists is there. I can't find it, using either python or IDL. But here is a shot...
;
; See my e-mails with Randy 18-19 Aug 2016.
;
; HBT 19-Aug-2016

common units

; Read the Level-2 Alice data

file_level2 =  '/Users/throop/Data/NH_Alice_Ring/O_RING_OC3/data/pluto/level2/ali/all/ali_0299391413_0x4b5_sci_1.fit'

im0 = mrdfits(file_level2,0,h1) ; PRIMARY
im1 = mrdfits(file_level2,1,h1) ; ERROR_IMAGE
im2 = mrdfits(file_level2,2,h2) ; WAVELENGTH_IMAGE  
im3 = mrdfits(file_level2,3,h3) ; PIXEL_LIST_TABLE
im4 = mrdfits(file_level2,4,h4) ; COUNT_RATE
im5 = mrdfits(file_level2,5,h5) ; HOUSEKEEPING_TABLE

print "# timesteps = " + st(sizex(im4))
print "# photons = " + st(sizex(im3))

a5 = mrdfits(file, 5)  

print a5

; Read the Level-1 Alice data

file_level1 = '/Users/throop/Data/NH_Alice_Ring/O_RING_OC3/data/pluto/level1/ali/all/ali_0299391413_0x4b5_eng_1.fit'

im0 = mrdfits(file_level1,0,h1) ; PRIMARY
im1 = mrdfits(file_level1,1,h1) ; PIXELLIST
im2 = mrdfits(file_level1,2,h2) ; COUNTRATE  
im3 = mrdfits(file_level1,3,h3) ; HOUSEKEEPING  ; Randy says that HK should be extension 2, but I think it is 3

count_rate = im3.count_rate + 2L^15

plot, im3.met+2L^31, im3.count_rate+2L^15, psym = 10  ; Randy says this should work

; Now try one more way of reading it. Interestingly, this finds a few structures that 
; are *not* found by mrdfits() or by python. But I'm not sure what they are. The length 
; is not really correlated with anything I'd expect.

rdfits_struct, file_level1, s3, exten=3

help, s3.tab3
help, s3.im0
help, s3.hdr0

