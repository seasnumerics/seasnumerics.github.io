// POV-Ray 3.6 / 3.7 Scene File "raster_plate_4.pov"
// author: Friedrich A. Lohmueller, Feb-2011
// email: Friedrich.Lohmueller_at_t-online.de
// homepage: http://www.f-lohmueller.de
//--------------------------------------------------------------------------
#version 3.6; // 3.7
global_settings { assumed_gamma 1.0 } 
#default{ finish{ ambient 0.1 diffuse 0.9 conserve_energy}}
//--------------------------------------------------------------------------
#include "colors.inc"
#include "textures.inc"
#include "glass.inc"
#include "metals.inc"
#include "golds.inc"
#include "stones.inc"
#include "woods.inc"
#include "shapes.inc"
#include "shapes2.inc"
#include "functions.inc"
#include "math.inc"
#include "transforms.inc"
//------------------------------------------------------------- Camera_Position, Camera_look_at, Camera_Angle
#declare Camera_Number = 0;


//--------------------------------------------------------------------------------------------------------<<<<
#switch ( Camera_Number )
#case (0)
  #declare Ultra_Wide_Angle_On = 0; 
  #declare Orthographic_On = 0; 
  #declare Camera_Position = < 1.50, 1.50,  1.50> ;  // front view                                                                 
  #declare Camera_Look_At  = < 0.00, 0.00,  0.00> ;
  #declare Camera_Angle    =  0 ;
  #declare Camera_Rotate = <0,0,0>;  
#break
#case (1)
  #declare Ultra_Wide_Angle_On = 0; 
  #declare Orthographic_On = 2; 
  #declare Camera_Position = < 1000,1700,250>;  // diagonal view < 1000,1300,250>
  #declare Camera_Look_At  = < 1000,250,250>;
  #declare Camera_Angle    =  65;
  #declare Camera_Rotate = <-10,-10,-25>;  
#break
#case (2)
  #declare Ultra_Wide_Angle_On = 0; 
  #declare Orthographic_On = 0; 
  #declare Camera_Position = < 0.00,10.00, -0.001> ;  // top view
  #declare Camera_Look_At  = < 0.00, 0.00,  0.000> ;
  #declare Camera_Angle    =  65 ;
  #declare Camera_Rotate = <0,0,0>;  
#break
#else
  #declare Ultra_Wide_Angle_On = 0; 
  #declare Orthographic_On = 0; 
  #declare Camera_Position = < 0.00, 1.00,-20.00> ;  // front view
  #declare Camera_Look_At  = < 0.00, 1.00,  0.00> ;
  #declare Camera_Angle    =  65 ;
  #declare Camera_Rotate = <0,0,0>;  
#break
#end // of "#switch ( Camera_Number )" -----------------------------
//-------------------------------------------------------------------------------------------------------<<<<

camera{ #if (Ultra_Wide_Angle_On = 1) ultra_wide_angle #end  
        #if (Orthographic_On = 1)     orthographic     #end  
        location Camera_Position
        right    x*image_width/image_height
        angle    Camera_Angle
        look_at  Camera_Look_At
        rotate   Camera_Rotate
      }
//------------------------------------------------------------------------------------------------------<<<<<

//------------------------------------------------------------------------
// sun -------------------------------------------------------------------
light_source{<2500,2500,2500> color White*0.9}           // sun light
light_source{ Camera_Position  color rgb<0.9,0.9,1>*0.1}  // flash light
background {color White}

//------------------------------ the Axes --------------------------------
//------------------------------------------------------------------------
#macro Axis_( AxisLen, Dark_Texture,Light_Texture)
 union{
    cylinder { <0,-AxisLen,0>,<0,AxisLen,0>,0.05
               texture{checker texture{Dark_Texture }
                               texture{Light_Texture}
                       translate<0.1,0,0.1>}
             }
    cone{<0,AxisLen,0>,0.2,<0,AxisLen+0.7,0>,0
          texture{Dark_Texture}
         }
     } // end of union
#end // of macro "Axis()"
//------------------------------------------------------------------------
#macro AxisXYZ( AxisLenX, AxisLenY, AxisLenZ, Tex_Dark, Tex_Light)
//--------------------- drawing of 3 Axes --------------------------------
#declare Text_Rotate = <20,-35,0>;
union{
#if (AxisLenX != 0)
 object { Axis_(AxisLenX, Tex_Dark, Tex_Light)   rotate< 0,0,-90>}// x-Axis
 text   { ttf "arial.ttf",  "x",  0.15,  0  texture{Tex_Dark}
          rotate Text_Rotate scale 0.5 translate <AxisLenX+0.15,0.3,-0.05> no_shadow }
#end // of #if
#if (AxisLenY != 0)
 object { Axis_(AxisLenY, Tex_Dark, Tex_Light)   rotate< 0,0,  0>}// y-Axis
 text   { ttf "arial.ttf",  "y",  0.15,  0  texture{Tex_Dark}
          rotate <Text_Rotate.x,0,0> scale 0.5 translate <-0.45,AxisLenY+0.20,-0.05> rotate <0,Text_Rotate.y,0> no_shadow }
#end // of #if
#if (AxisLenZ != 0)
 object { Axis_(AxisLenZ, Tex_Dark, Tex_Light)   rotate<90,0,  0>}// z-Axis
 text   { ttf "arial.ttf",  "z",  0.15,  0  texture{Tex_Dark}
          rotate Text_Rotate scale 0.65 translate <-0.75,0.1,AxisLenZ+0.10> no_shadow }
#end // of #if
} // end of union
#end// of macro "AxisXYZ( ... )"

//------------------------------------------------------------------------

#declare Texture_A_Dark  = texture {
                               pigment{ color rgb<1,0.45,0>}
                               finish { phong 1}
                             }
#declare Texture_A_Light = texture {
                               pigment{ color rgb<1,1,1>}
                               finish { phong 1}
                             }

//object{ AxisXYZ( 4.5, 4.5, 2.7, Texture_A_Dark, Texture_A_Light) scale 0.5}
//-------------------------------------------------- end of coordinate axes

//---------------------------- objects in scene ----------------------------
camera { 
direction z 
location <0, 6, -15> 
look_at <0, 2, 0> 
} 
//union{
#include "seafloor2.pov"
//}
//union{
#include "Xmemb.inc"
object {m_dat texture{ pigment{ color rgbt<0.0, 0.9, 0.9, 0.5>}}}
//#include "shapes3.inc"
//}
