# Copyright  2014 C.H. Wu <chenhung@cims.nyu.edu>. All rights reserved.
# Date: 2014-08-31 02:09:40
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY AUTHOR AND CONTRIBUTORS ``AS IS'' AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED.  IN NO EVENT SHALL AUTHOR OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.

"""Program for converting an vtk file into a POV-ray mesh or mesh2."""

import argparse
import sys
import time
import pyvtk
import utils
from os.path import basename

__version__ = '$Revision: 1.0 $'[11:-2]

def vtk_to_mesh1(name, vertices, triangles):
    """Creates a POV-ray mesh description from facet data.

    :name: The name of the object.
    :vertices: An (N,3) array containing the vertex data.
    :triangles: An (N,3) array containing the triangle data 
    :returns: a string representation of a POV-ray mesh object.
    """
    lines = ["# declare m_{} = mesh {{".format(name.replace(' ', '_'))]
    sot = "  triangle {"
    # The indices sequence 1, 0, 2 is used because of the difference between
    # the VTK coordinate system and that used in POV-ray.
    fc = "    <{0}, {1}, {2}>,"
    for (a, b, c) in triangles:
        lines += [sot, fc.format(*vertices[a]), fc.format(*vertices[b]),
                  fc.format(*vertices[c])[:-1], "  }"]
    lines += ['}']
    return '\n'.join(lines)

def vtk_to_sphere(name, particles, radius):
    """Creates a POV-ray mesh description from facet data.

    :name: The name of the object.
    :vertices: An (N,3) array containing the vertex data.
    :triangles: An (N,3) array containing the triangle data 
    :returns: a string representation of a POV-ray mesh object.
    """
    lines = ["# declare m_{} = union {{".format(name.replace(' ', '_'))]
    sot = "  sphere {"
    # The indices sequence 1, 0, 2 is used because of the difference between
    # the VTK coordinate system and that used in POV-ray.
    fc = "    <{0}, {1}, {2}>,"
    for prt in particles:
        lines += [sot, fc.format(*prt), radius, "  }"]
    lines += ['}']
    return '\n'.join(lines)


def mesh1(name, vertices):
    """Creates a POV-ray mesh description from facet data.

    :name: The name of the object.
    :vertices: An (N,3) numpy array containing the vertex data.
    :returns: a string representation of a POV-ray mesh object.
    """
    facets = vertices.reshape((-1, 3, 3))
    lines = ["# declare m_{} = mesh {{".format(name.replace(' ', '_'))]
    sot = "  triangle {"
    # The indices sequence 1, 0, 2 is used because of the difference between
    # the STL coordinate system and that used in POV-ray.
    fc = "    <{1}, {0}, {2}>,"
    for (a, b, c) in facets:
        lines += [sot, fc.format(*a), fc.format(*b),
                  fc.format(*c)[:-1], "  }"]
    lines += ['}']
    return '\n'.join(lines)


def read_solutes(filename):
    """Reads a solute vtk file, returns the points and the name.
    The POINT_DATA information is discarded since it is not goanna use.

    :name: path of the vtk file to read
    :returns: a array of the shape (N, 3) containing the points and the 
    name of the object as given in the vtk file.
    """
    name = basename(filename)[:-8] # drop the .vtk extension
    lines = [line.strip() for line in open(filename)]
    num_solutes = int(lines[5].split()[1])
    particles = []
    for l in range(6, 6+num_solutes):
        particles.append(tuple([float(p) for p in lines[l].split()]))
    return particles, name

def read_memb(filename):
    """Reads a vtk file, returns the vertices, triangles, and the name.
    The CELL_TYPES information is discarded since it is not goanna use.

    :name: path of the vtk file to read
    :returns: a numpy array of the shape (N, 3) containing the vertices
    of the facets, and the name of the object as given in the vtk file.
    """
    vtk_data = pyvtk.VtkData(filename)
    vertices = vtk_data.structure.points
    triangles = vtk_data.structure.triangle
    name = basename(filename)[:-8] # drop the 0005.vtk extension 
    if vertices is None:
        raise ValueError('not a valid VTK file.')
    return vertices, triangles, name

def main(argv):
    """Main program.

    :argv: command line arguments (without program name!)
    """
    msg = utils.Msg()
    parser = argparse.ArgumentParser(description=__doc__)
    argtxt = 'generate a mesh2 object (slow on big files)'
    parser.add_argument('-2,' '--mesh2', action='store_true',
                        help=argtxt, dest='mesh2')
    parser.add_argument('-v', '--version', action='version',
                        version=__version__)
    parser.add_argument('file', nargs='*', help='one or more file names')
    args = parser.parse_args(argv)
    if not args.file:
        parser.print_help()
        sys.exit(0)
    for fn in args.file:
        bfn = basename(fn)
        msg.say('Starting file "{}"'.format(fn))
        if not fn.lower().endswith('.vtk'):
            w = 'The file "{}" is probably not a vtk file, skipping.'
            msg.say(w.format(fn))
            continue
        try:
            msg.say('Reading vtk')
            #if bfn.lower().startswith('solute'):
            #    particles, name = read_solutes(fn)
            #if bfn.lower().startswith('xmemb'):
            vertices, triangles, name = read_memb(fn)
            outfn = utils.outname(fn, '.inc')
        except Exception as e:  # pylint: disable=W0703
            msg.say('Error;', e)
            continue
        outs = "// Generated by vtk2pov {}\n".format(__version__)
        outs = "// on {}.\n".format(time.asctime())
        outs += "// Source file name: '{}'\n".format(fn)
        if args.mesh2:
            msg.say('Generating mesh2 data')
            outs += mesh2(name, vertices)
        else:
            msg.say('Generating mesh data')
            if bfn.lower().startswith('solute'):
                outs += vtk_to_sphere(name, particles, '7')
            #if bfn.lower().startswith('xmemb'): 
            outs += vtk_to_mesh1(name, vertices, triangles)
        try:
            with open(outfn, 'w+') as of:
                msg.say('Writing output file "{}"'.format(outfn))
                of.write(outs)
        except:
            msg.say('Cannot write output file "{}"'.format(outfn))
            continue
        msg.say('Done with file "{}"'.format(fn))


if __name__ == '__main__':
    main(sys.argv[1:])

