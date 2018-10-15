#!/usr/bin/perl

# Total output frames
$maxf=2400;

# Options controlling the POV-Ray rendering quality
$pov_opts="+H768 +W1024 +A0.001 +R9 -J";

# Load in the solution that was previously stored using chua.py. Only store the rotated solutions
open A,"chua_sol.dat";
$n=0;
while(<A>) {
    @a=split;
    $x[$n]=$a[4]*1.1;
    $y[$n]=$a[5]*1.1;
    $z[$n++]=$a[6];
}
close A;
$n--;

mkdir "chua_frames";
$t=nonlin($maxf);
foreach $s (0..$maxf) {

    # Compute how far of the solution trace to save.
    $u=$n*nonlin($s)/$t;
    $iu=int $u;
    $fu=$u-$iu;

    # Compute the rotation angle
    $rot=360*$s/1200.;

    # Open the POV-Ray template and replace the markers like "AAA", etc.
    # with the corresponding parameter
    open A,"c_template.pov" or die "Can't open POV template\n";
    open B,">render.pov" or die "Can't open temporary file\n";
    while(<A>) {
        if(/CHUA/) {

            # Print a sphere at the very start of the trace
            printf B "sphere{<%.7g,%.7g,%.7g>,r %s finish{f0}}\n",$x[0],$y[0],$z[0],col(0);

            # Loop over the tracer segments
            foreach $i (0..$iu) {
                if($i==$iu) {

                    # The final segment may be a partial segment. Compute its
                    # end point using linear interpolation. If the segment is
                    # very small, then just ignore it.
                    break if $fu<1e-5;
                    $xx=(1-$fu)*$x[$i]+$fu*$x[$i+1];
                    $yy=(1-$fu)*$y[$i]+$fu*$y[$i+1];
                    $zz=(1-$fu)*$z[$i]+$fu*$z[$i+1];
                    $ha=0.5*$fu;
                } else {

                    # Otherwise, print the complete segment
                    $xx=$x[$i+1];
                    $yy=$y[$i+1];
                    $zz=$z[$i+1];
                    $ha=0.5;
                }

                # POV-Ray gives an error if cylinders have the same start and
                # end point. Skip this cylinder if the two positions match.
                $pos1=sprintf "<%.7g,%.7g,%.7g>",$x[$i],$y[$i],$z[$i];
                $pos2=sprintf "<%.7g,%.7g,%.7g>",$xx,$yy,$zz;
                printf B "cylinder{$pos1,$pos2,r %s finish{f0}}\n",col(($i+$ha)/$n) unless $pos1 eq $pos2;

                # Print the sphere at the end point of this segment
                printf B "sphere{$pos2,r %s finish{f0}}\n",col($i/$n);
            }
            next;
        }
        s/ROT/$rot/;
        print B;
    }
    close A;
    close B;

    # System call to run POV-Ray and output a frame of the movie
    print "Frame $s: $u\n";
    $of=sprintf "chua_frames/fr_%04d.png",$s;
    system "nice -n 19 povray $pov_opts +O$of render.pov >/dev/null 2>/dev/null";
}

# Nonlinear function that governs how fast the tracer is plotted. It is
# designed to speed up as the movie progresses.
sub nonlin {
    return @_[0]*(1+sqrt(@_[0]));
}

# Nonlinear scaling to apply to the color palette
sub rfunc {
    return @_[0]*@_[0]*(3-2*@_[0]);
}

# Subroutine to return a rainbow color based on an input from 0 to 1.
sub col {
    $ang=@_[0]*6.28318530717959;
    $re=0.5+0.5*cos($ang);
    $gr=0.5+0.5*cos($ang-2.09439510239329549230842892219);
    $bl=0.5+0.5*cos($ang-4.18879020478649098461685784436);

    # Return the POV-Ray pigment command. Since each color channel only has 256
    # entries, there is little need to store more than four decimal places of
    # accuracy.
    return sprintf "pigment{rgb <%.4f,%.4f,%.4f>}",rfunc($re),rfunc($gr),rfunc($bl);
}
