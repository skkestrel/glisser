# -*- coding: UTF-8 -*-
import numpy as np

# for Mercury and Swift
mu_swift = 4*np.pi**2
mu_mercury = mu_swift

# for GLISSE
mu_glisse = 4*np.pi**2 / (365.25)**2

m_sun = 1.98847e30
year = 365.25 * 24 * 3600
au = 149597870700
vel_scale = 2*np.pi/year*au


def e2c(inp, mu_sun):
    """
    input: {{a, e, I, g, n, f}} \\
    output: {{x, y, z, u, v, w}}(inerital)

    a = semi-major axis(in AU) \\
    e = eccentricity \\
    I = inclination (degrees) \\
    g = argument of pericentre (degrees) or omega \\
    n = longitude of the ascending node (degrees) or Omega \\
    f = true anomaly (degrees)
    """
    a, e, I, g, n, f = inp[0], inp[1], np.deg2rad(inp[2]), np.deg2rad(
        inp[3]), np.deg2rad(inp[4]), np.deg2rad(inp[5])

    p = a * np.abs(1-e**2)

    sini, cosi, sing, cosg, sinn, cosn = np.sin(I),  np.cos(
        I), np.sin(n), np.cos(n), np.sin(g), np.cos(g)

    HVector = np.array([sini * sing, -sini * cosg, cosi])
    PVector = np.array([cosg * cosn - sing * sinn * cosi,
                        sing * cosn + cosg * sinn * cosi,
                        sinn * sini])
    QVector = np.cross(HVector, PVector)

    r = p / (1.0 + e * np.cos(f))
    out = [0, 0, 0, 0, 0, 0]
    for i in range(3):
        out[i] = r * (np.cos(f) * PVector[i] + np.sin(f) * QVector[i])
        out[3 + i] = np.sqrt(mu_sun / p) * \
            (-np.sin(f) * PVector[i] + (np.cos(f) + e) * QVector[i])
    return out


def c2e(inp, mu_sun):
    """
    input: {{x,y,z,u,v,w}} (inerital) \\
    output: {{a,e,I,g,n,f}}

    a = semi-major axis (in AU) \\
    e = eccentricity \\
    I = inclination (degrees) \\
    g = argument of pericentre (degrees) or omega \\
    n = longitude of the ascending node (degrees) or Omega \\
    f = true anomaly (degrees) \\
    """
    R, V = np.array([inp[0], inp[1], inp[2]]), np.array(
        [inp[3], inp[4], inp[5]])
    radius, vel = np.linalg.norm(R), np.linalg.norm(V)
    vr = np.dot(R, V) / radius
    unitR = R/radius
    unitV = V/vel
    hvector = np.cross(R, V)
    hnorm = np.linalg.norm(hvector)
    unith = hvector/hnorm
    temp1 = vel**2 - mu_sun / radius
    temp2 = radius * vr
    evector = [0, 0, 0]
    for i in range(3):
        evector[i] = (temp1 * R[i] - temp2 * V[i]) / mu_sun
    e = np.linalg.norm(evector)
    isCircle = (np.abs(e) <= 1e-15)
    a = hnorm**2 / (mu_sun * (1 - e**2))
    I = np.arccos(unith[2])
    unitN = np.array([-unith[1], unith[0], 0])
    if (np.linalg.norm(unitN) == 0):
        n = 0
        if (isCircle):
            g = 0
            f = np.arctan2(unitR[1] * unith[2], unitR[0])
        else:
            unite = evector/e
            temp = np.cross(unite, unitR)
            g = np.arctan2(unite[1] * unith[2], unite[0])
            f = np.arctan2(np.dot(unith, temp), np.dot(unite, unitR))
    else:
        temp = np.cross(unitN, unitR)
        n = np.arctan2(unith[0], -unith[1])
        f = np.arctan2(np.dot(unith, temp), np.dot(unitN, unitR))
        if (isCircle):
            g = 0
        else:
            unite = evector/e
            temp = np.cross(unitN, unite)
            g = np.arctan2(np.dot(unith, temp), np.dot(unite, unitN))
            f = f - g

    if (g < 0):
        g += 2*np.pi
    if (n < 0):
        n += 2*np.pi
    if (f < 0):
        f += 2*np.pi

    I, g, n, f = np.rad2deg(I), np.rad2deg(g), np.rad2deg(n), np.rad2deg(f)
    return [a, e, I, g, n, f]


def true2mean(theta, e):
    # input: radians
    # output: radians
    if (theta == 0.0):
        return 0.0
    E = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(theta / 2))
    result = (E - e * np.sin(E))
    if result < 0:
        return result + 2*np.pi
    else:
        return result


def mean2true(M, e):
    # input: radians
    # output: radians
    acc = 1e-14
    if (M == 0.0):
        return 0.0
    if (M < np.pi and M > 0):
        E = M + e / 2
    else:
        E = M - e / 2
    incr = 1
    count = 0
    while (np.abs(incr) > acc):
        incr = keplerIteration(E, e, M)
        E += incr
    result = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))
    return result


def keplerIteration(E, e, M):
    return -(E - e * np.sin(E) - M) / (1 - e * np.cos(E))


def scale2real(ele):
    return [ele[0]*au, ele[1]*au, ele[2]*au, ele[3]*vel_scale, ele[4]*vel_scale, ele[5]*vel_scale]


def scale2non(ele):
    return [ele[0]/au, ele[1]/au, ele[2]/au, ele[3]/vel_scale, ele[4]/vel_scale, ele[5]/vel_scale]
