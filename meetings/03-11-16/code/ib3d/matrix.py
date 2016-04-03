# file: matrix.py
# vim:fileencoding=utf-8
#
# Copyright © 2013,2014 R.F. Smith <rsmith@xs4all.nl>. All rights reserved.
# Created: 2013-07-28 02:07:00 +0200
# $Date: 2014-06-15 17:12:40 +0200 $
# $Revision: 3.2 $
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
# THIS SOFTWARE IS PROVIDED BY AUTHOR AND CONTRIBUTORS "AS IS" AND
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

"""3D homogeneous coordinates matrix functions in a
right-handed coordinate system."""

import math
import numpy as np

__version__ = '$Revision: 3.2 $'[11:-2]


def I():
    """Returns the 4x4 identity matrix."""
    return np.identity(4, np.float32)


def trans(vec):
    """Returns a 4x4 homogeneous coordinates translation matrix along vec.

    :vec: 3D translation vector
    """
    rv = I()
    rv[0, 3] = vec[0]
    rv[1, 3] = vec[1]
    rv[2, 3] = vec[2]
    return rv


def mul(begin, *rest):
    """Returns the multiplication of the 4x4 matrix arguments.

    :begin, rest: 4x4 matrices
    """
    rv = np.copy(begin)
    for r in rest:
        rv = np.dot(rv, r)
    return rv


def concat(*args):
    """Concatenate the transforms. This is actually a multiplication
    of the arguments in reversed order.

    :*args: 4x4 matrices
    :returns: a combined matrix
    """
    rv = np.copy(args[-1])
    rest = list(reversed(args[:-1]))
    for r in rest:
        rv = np.dot(rv, r)
    return rv


def rotx(angle):
    """Returns the 4x4 homogeneous coordinates matrix for rotation around the
    X axis.

    :angle: rotation angle in degrees
    """
    rad = math.radians(angle)
    c = math.cos(rad)
    s = math.sin(rad)
    return np.array([[1.0, 0.0, 0.0, 0.0],
                     [0.0,   c,  -s, 0.0],
                     [0.0,   s,   c, 0.0],
                     [0.0, 0.0, 0.0, 1.0]], np.float32)


def roty(ang):
    """Returns the 4x4 homogeneous coordinates matrix for rotation around the
    Y axis.

    :angle: rotation angle in degrees
    """
    rad = math.radians(ang)
    c = math.cos(rad)
    s = math.sin(rad)
    return np.array([[c, 0.0,   s, 0.0],
                     [0.0, 1.0, 0.0, 0.0],
                     [-s, 0.0,   c, 0.0],
                     [0.0, 0.0, 0.0, 1.0]], np.float32)


def rotz(ang):
    """Returns the 4x4 homogeneous coordinates matrix for rotation around the
    Z axis.

    :angle: rotation angle in degrees
    """
    rad = math.radians(ang)
    c = math.cos(rad)
    s = math.sin(rad)
    return np.array([[c,  -s, 0.0, 0.0],
                     [s,   c, 0.0, 0.0],
                     [0.0, 0.0, 1.0, 0.0],
                     [0.0, 0.0, 0.0, 1.0]]), np.float32


def rot(axis, angle):
    """Returns the 4x4 homogeneous coordinates matrix for rotation around the
    an arbitrary axis by an arbitrary angle.

    :axis: rotation axis
    :angle: rotation angle in degrees
    """
    ax = np.require(axis[0:3], np.float32)
    l = np.linalg.norm(ax)
    if l == 0.0:
        raise ValueError('axis cannot have length 0')
    elif not l == 1.0:
        ax /= l
    ux, uy, uz = ax
    a = math.radians(angle)
    c = math.cos(a)
    s = math.sin(a)
    uc = np.array([[0, -uz, uy],
                   [uz, 0, -ux],
                   [-uy, ux, 0]], np.float32)
    ut = np.array([[ux*ux, ux*uy, ux*uz],
                   [ux*uy, uy*uy, uy*uz],
                   [ux*uz, uy*uz, uz*uz]], np.float32)
    m = np.identity(3, np.float32)*c + uc*s + ut*(1.0 - c)
    m = np.vstack((m, np.array([0, 0, 0], np.float32)))
    m = np.column_stack((m, np.array([0, 0, 0, 1], np.float32)))
    return m


def scale(x=1, y=1, z=1):
    """Returns a scaling matrix

    :x: scale factor for x direction
    :y: scale factor for y direction
    :z: scale factor for z direction
    """
    rv = I()
    rv[0, 0], rv[1, 1], rv[2, 2] = float(x), float(y), float(z)
    return rv


def lookat(eye, center, up):
    """Create a viewing matrix

    :eye: 3D point where the viewer is located
    :center: 3D point that the eye looks at
    :up: 3D upward direction
    :returns: view matrix
    """
    eye = np.array(eye, np.float32)
    center = np.array(center, np.float32)
    up = np.array(up, np.float32)
    F = center - eye
    f = F/np.linalg.norm(F)
    S = np.cross(f, up)
    s = S/np.linalg.norm(S)
    u = np.cross(s, f)
    rv = np.array([[s[0],  s[1],  s[2], -eye[0]],
                   [u[0],  u[1],  u[2], -eye[1]],
                   [-f[0], -f[1], -f[2], -eye[2]],
                   [0,     0,     0,       1]], np.float32)
    return rv


def ortho(xyscale):
    """Creates a simple orthographic projection matrix.
    :xyscale: scaling factor for x and y
    :returns: orthographic projection matrix
    """
    rv = I()
    rv[0, 0], rv[1, 1] = xyscale, xyscale
    rv[2, 2] = 0
    return rv


def perspective(fovy, width, height, near, far):
    """Create a prespective projection matrix.

    :fovy: field of view in y direction, in degrees
    :width: width of the viewport
    :height: height of the viewport
    :znear: near clipping plane (> 0)
    :zfar: far clipping plane (> 0)
    :returns: projection matrix
    """
    aspect = float(width)/float(height)
    f = 1/math.tan(math.radians(float(fovy))/2)
    near = float(near)
    far = float(far)
    d = near - far
    rv = np.array([[f/aspect, 0, 0, 0],
                   [0, f, 0, 0],
                   [0, 0, (far+near)/d, 2*far*near/d],
                   [0, 0, -1, 0]], np.float32)
    return rv
