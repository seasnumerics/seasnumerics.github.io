// Persistence Of Vision raytracer version 3.5 sample file. 
// Caustics example 

global_settings { assumed_gamma 2.2 } 

#include "colors.inc" 
#include "textures.inc" 

light_source { <0, 50, 0> color White } 

camera { 
  direction z 
  location <1.5, 1.5, 1.5> 
  look_at <0, 0, 0> 
} 

// The sea floor 
plane { y, 0 
  pigment { Gray60 } 
  finish { ambient 0.1 diffuse 0.7 } 
} 

// The water surface 
plane { y, 10 
  hollow on 
  pigment { red 0.7 green 0.7 blue 1.0 filter 0.9 } 
  finish {reflection 0.7 } 
  interior { ior 1.1 caustics 1.0 } 
  translate <5, 0, -10> 
  normal { bumps 0.5 } 
}

union{
  #include "Xmemb.inc"
  object {m_ material{texture{ pigment{ color rgbt<1,1,1,0.5>}
        normal { bumps 0.5 scale 0.025 }
        finish { phong 1.0 }
      } // end of texture
      interior{ ior 1.5
        caustics 0.25
      } // end of interior
    } 
  }
}

