Definition of '2014 MU69 Sunflower' reference frames.

Developed for New Horizons KEM MU69 flyby & hazard search.

The frames defined here follow the system proposed in Mark Showalter's
4-Jan-2018 'Coordinate Frames Strawman Proposal,' as implemented by Henry
Throop

Frames defined in this file are:

  '2014_MU69_SUNFLOWER_ROT'  
  '2014_MU69_SUNFLOWER_INERT'
  '2014_MU69_TUNACAN_ROT'  
  '2014_MU69_TUNACAN_INERT'
  '2014_MU69_ORT4_1'

Sunflower Frames
---
The first is a *non-inertial* coordinate frame such that:

  * +Y_SUNFLOWER points from MU69 center, away from the (apparent) Sun center;
  * +Z_SUNFLOWER points in the direction of MU69's orbital pole.
  * +X_SUNFLOWER completes the right-handed system.

This is a non-inertial frame. It rotates once every time MU69 
orbits the Sun (~296 years).

We can also define a related *inertial* coordinate frame, which is identical
except:

  * +Y_SUNFLOWER points from MU69 center, to the fixed RA/Dec in the anti-sun
    direction at 5-Dec-2018 00:00:00.

These are both Type-5 frames, aka a 'Parameterized Dynamic Frame.'

Primary vector: Define -Y to point in the MU69-Sun direction.

Secondary vector: Defined to the orbital pole of MU69. This is taken from value
calculated by Simon Porter 15-Jan-2018. The orbit pole position may be refined
during the MU69 approach period, but any changes will be trivially small.

    MU69 osculating orbit pole on 2019-01-01T00:00:00
    ICRF RA: 272.426110231801, Dec: 68.831520928192

Tunacan Frames
---

These frames are exactly the same as the Sunflower frames, but rotated 90 deg 
s.t. the XZ orbital plane around MU69 is roughly in the Solar System's
equatorial plane, rather than facing the Sun. Rings in this configuration would
look more like a 'tuna can' than a 'sunflower.'

The tuna can frame is incorporated into the backplane generation for NH Hazard
work. I don't expect to use it very broadly.

NB: The reference to 'J2000' in FRAME_*_RELATIVE specifies a computation frame
for SPICE's internal usage, but the results returned are nearly identical
whether this is J2000, B1950, or any other inertial reference frame.

ORT4 Frames
---

This is a one-off frame which has the +Y axis pointing toward a specific
RA/Dec, so as to match the obsevations. The position is determined from ORT4
data, using NH_ORT_FIND_RING_POLE.PY.

The second axis determines the zero-point of the longitude system, and is
basically arbitrary. I set it to (0,0).

ORT4_1: Pole points at (275, -56) deg.

ORT4_n: Reserved for future use, and different pole positions.

---

Henry Throop, 16-Jan-2018. Initial release to Hazard team. Please let me know
                           any issues.

Henry Throop, 25-May-2018. Added MU69 Tuna Can planes.

Henry Throop, 22-Jun-2018. Added MU69 ORT4_1 frame.

---
Define ROTATING (NON-INERTIAL) reference frame.
---

\begindata

  FRAME_2014_MU69_SUNFLOWER_ROT		= 1234580
  FRAME_1234580_NAME			= '2014_MU69_SUNFLOWER_ROT'
  FRAME_1234580_CLASS			= 5
  FRAME_1234580_CENTER			= '2014 MU69'
  FRAME_1234580_CLASS_ID		= 1234580
  FRAME_1234580_RELATIVE		= 'J2000'
  FRAME_1234580_DEF_STYLE		= 'PARAMETERIZED'
  FRAME_1234580_FAMILY			= 'TWO-VECTOR'
  FRAME_1234580_ROTATION_STATE		= 'ROTATING'

  FRAME_1234580_PRI_AXIS		= '-Y'
  FRAME_1234580_PRI_VECTOR_DEF		= 'OBSERVER_TARGET_POSITION'
  FRAME_1234580_PRI_OBSERVER		= '2014 MU69'
  FRAME_1234580_PRI_TARGET		= 'SUN'
  FRAME_1234580_PRI_ABCORR		= 'LT'
  FRAME_1234580_PRI_FRAME		= 'J2000'

  FRAME_1234580_SEC_AXIS		= '+Z'
  FRAME_1234580_SEC_VECTOR_DEF          = 'CONSTANT'
  FRAME_1234580_SEC_SPEC                = 'RA/DEC'
  FRAME_1234580_SEC_UNITS               = 'DEGREES'
  FRAME_1234580_SEC_RA                  = 272.426110231801
  FRAME_1234580_SEC_DEC                 = 68.831520928192
  FRAME_1234580_SEC_FRAME		= 'J2000'

\begintext

---
Define INERTIAL reference frame.
---

\begindata

  FRAME_2014_MU69_SUNFLOWER_INERT	= 1234581
  FRAME_1234581_NAME                    = '2014_MU69_SUNFLOWER_INERT'
  FRAME_1234581_CLASS                   = 5
  FRAME_1234581_CENTER                  = '2014 MU69'
  FRAME_1234581_CLASS_ID                = 1234581
  FRAME_1234581_RELATIVE                = 'J2000'
  FRAME_1234581_DEF_STYLE               = 'PARAMETERIZED'
  FRAME_1234581_FAMILY                  = 'TWO-VECTOR'
  FRAME_1234581_FREEZE_EPOCH            = @5-DEC-2018/00:00:00

  FRAME_1234581_PRI_AXIS                = '-Y'
  FRAME_1234581_PRI_VECTOR_DEF          = 'OBSERVER_TARGET_POSITION'
  FRAME_1234581_PRI_OBSERVER            = '2014 MU69'
  FRAME_1234581_PRI_TARGET              = 'SUN'
  FRAME_1234581_PRI_ABCORR              = 'LT'
  FRAME_1234581_PRI_FRAME		= 'J2000'
 
  FRAME_1234581_SEC_AXIS                = '+Z'
  FRAME_1234581_SEC_VECTOR_DEF          = 'CONSTANT'
  FRAME_1234581_SEC_SPEC                = 'RA/DEC'
  FRAME_1234581_SEC_UNITS               = 'DEGREES'
  FRAME_1234581_SEC_RA                  = 272.426110231801
  FRAME_1234581_SEC_DEC                 = 68.831520928192
  FRAME_1234581_SEC_FRAME		= 'J2000'

\begintext

** MU69 TUNACAN FRAME **

Now we define the 'Tuna Can' frame. It is called this because a ring in the XZ
plane of this body will resemble a tunacan

I am defining an entirely new frame, s.t. my backplane code (which assumes a
ring in the XZ plane) can generate backplanes for a Tuncan orbit as easily as
for a Sunflower orbit.

This is coordinate system that is rotated 90 degrees toward the sun. It is an
edge-on sunflower orbit. This orbit and kernel will be used primarily for
generating a backplane for the ORT. I don't anticipate using the tunacan 
frame for anything else -- e.g., not in any of the ORT steps, or for DPH /
Mehoke calculations, etc.

  * +Z_TUNACAN   points from MU69 center, away from the (apparent) Sun center;
  * +Y_TUNACAN   points in the direction of MU69's orbital pole.
  * +X_TUNACAN   completes the right-handed system.

---
Define ROTATING (NON-INERTIAL) reference frame, MU69_TUNACAN
---

\begindata

  FRAME_2014_MU69_TUNACAN_ROT		= 1234582
  FRAME_1234582_NAME			= '2014_MU69_TUNACAN_ROT'
  FRAME_1234582_CLASS			= 5
  FRAME_1234582_CENTER			= '2014 MU69'
  FRAME_1234582_CLASS_ID		= 1234582
  FRAME_1234582_RELATIVE		= 'J2000'
  FRAME_1234582_DEF_STYLE		= 'PARAMETERIZED'
  FRAME_1234582_FAMILY			= 'TWO-VECTOR'
  FRAME_1234582_ROTATION_STATE		= 'ROTATING'

  FRAME_1234582_PRI_AXIS		= '-Z'
  FRAME_1234582_PRI_VECTOR_DEF		= 'OBSERVER_TARGET_POSITION'
  FRAME_1234582_PRI_OBSERVER		= '2014 MU69'
  FRAME_1234582_PRI_TARGET		= 'SUN'
  FRAME_1234582_PRI_ABCORR		= 'LT'
  FRAME_1234582_PRI_FRAME		= 'J2000'

  FRAME_1234582_SEC_AXIS		= '+Y'
  FRAME_1234582_SEC_VECTOR_DEF          = 'CONSTANT'
  FRAME_1234582_SEC_SPEC                = 'RA/DEC'
  FRAME_1234582_SEC_UNITS               = 'DEGREES'
  FRAME_1234582_SEC_RA                  = 272.426110231801
  FRAME_1234582_SEC_DEC                 = 68.831520928192
  FRAME_1234582_SEC_FRAME		= 'J2000'

\begintext

---
Define INERTIAL reference frame, MU69_TUNACAN
---

\begindata

  FRAME_2014_MU69_TUNACAN_INERT	        = 1234583
  FRAME_1234583_NAME                    = '2014_MU69_TUNACAN_INERT'
  FRAME_1234583_CLASS                   = 5
  FRAME_1234583_CENTER                  = '2014 MU69'
  FRAME_1234583_CLASS_ID                = 1234583
  FRAME_1234583_RELATIVE                = 'J2000'
  FRAME_1234583_DEF_STYLE               = 'PARAMETERIZED'
  FRAME_1234583_FAMILY                  = 'TWO-VECTOR'
  FRAME_1234583_FREEZE_EPOCH            = @5-DEC-2018/00:00:00

  FRAME_1234583_PRI_AXIS                = '-Z'
  FRAME_1234583_PRI_VECTOR_DEF          = 'OBSERVER_TARGET_POSITION'
  FRAME_1234583_PRI_OBSERVER            = '2014 MU69'
  FRAME_1234583_PRI_TARGET              = 'SUN'
  FRAME_1234583_PRI_ABCORR              = 'LT'
  FRAME_1234583_PRI_FRAME		= 'J2000'
 
  FRAME_1234583_SEC_AXIS                = '+Y'
  FRAME_1234583_SEC_VECTOR_DEF          = 'CONSTANT'
  FRAME_1234583_SEC_SPEC                = 'RA/DEC'
  FRAME_1234583_SEC_UNITS               = 'DEGREES'
  FRAME_1234583_SEC_RA                  = 272.426110231801
  FRAME_1234583_SEC_DEC                 = 68.831520928192
  FRAME_1234583_SEC_FRAME		= 'J2000'

\begintext

---
Define 2014_MU69_ORT4_1 frame. This is just a fixed frame that points to a
given RA/Dec.
---

\begindata

  FRAME_2014_MU69_ORT4_1	        = 1234584
  FRAME_1234584_NAME                    = '2014_MU69_ORT4_1'
  FRAME_1234584_CLASS                   = 5
  FRAME_1234584_CENTER                  = '2014 MU69'
  FRAME_1234584_CLASS_ID                = 1234584
  FRAME_1234584_RELATIVE                = 'J2000'
  FRAME_1234584_DEF_STYLE               = 'PARAMETERIZED'
  FRAME_1234584_FAMILY                  = 'TWO-VECTOR'
  FRAME_1234584_FREEZE_EPOCH            = @5-DEC-2018/00:00:00

  FRAME_1234584_PRI_AXIS                = '-Y'
  FRAME_1234584_PRI_VECTOR_DEF          = 'CONSTANT'
  FRAME_1234584_PRI_SPEC                = 'RA/DEC'
  FRAME_1234584_PRI_UNITS               = 'DEGREES'
  FRAME_1234584_PRI_ABCORR              = 'LT'
  FRAME_1234584_PRI_FRAME		= 'J2000'
  FRAME_1234584_PRI_RA                  = 275
  FRAME_1234584_PRI_DEC                 = -56
 
  FRAME_1234584_SEC_AXIS                = '+Z'
  FRAME_1234584_SEC_VECTOR_DEF          = 'CONSTANT'
  FRAME_1234584_SEC_SPEC                = 'RA/DEC'
  FRAME_1234584_SEC_UNITS               = 'DEGREES'
  FRAME_1234584_SEC_RA                  = 000
  FRAME_1234584_SEC_DEC                 = 00 
  FRAME_1234584_SEC_FRAME		= 'J2000'
\begintext
