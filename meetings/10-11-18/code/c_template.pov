#version 3.6;

// Specify a right-handed coordinate system in which the z-axis points upwards
camera {
    orthographic
	location <0,-50,0>
	sky z
	right -4.5*image_width/image_height*x
	up 4.5*z
	look_at <0,0,0>
}

// Two lights with slightly different colors
light_source{<-20,-30,30> color rgb <0.77,0.75,0.75>}
light_source{<25,-32,12> color rgb <0.38,0.40,0.40>}

#declare r=0.006;
#declare f0=finish{reflection 0.02 specular 0.35 ambient 0.45}

union {
CHUA
	rotate <0,0,ROT>
}
