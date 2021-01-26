GLISSER - GPU Long-term Integrator for Solar System Evolution (Regularized)
===========================================================================

Introduction
------------

GLISSER is an integrator which makes use of a GPU card to perform long-term and long-scale integrations of massless test particles in a solar system.

### Functions

GLISSER is capable of:

*   GPU-only integration using either an interpolated planetary history file or planetary initial conditions, which does not support particle close encounters
*   GPU-CPU hybrid integration using an interpolated planetary history file, which supports particle close encounters

### Requirements

*   A CUDA-enabled GPU (CUDA 8 recommended)
*   Linux

Installation
------------

### Compiling the main code

GLISSER functions as a standalone code. Simply type `make` in the project root, which will generate the executable `bin/glisse`. It may be necessary to edit the variable `CUDA_ARCH` in `Makefile` to match your system’s GPU architecture.

### Compiling the addon SWIFT code

Navigate into `swift/` and run `./CompileAll.sh`, which generates the executable `swift/main/swift_glisse_ext`. It may be necessary to edit the variables `NPLMAX` and `NTPMAX` in `swift/main/swift.inc` to increase the maximum number of particles or planets that SWIFT can integrate at once.

> [SWIFT](https://www.boulder.swri.edu/~hal/swift.html) is developed by Harold Levison and Martin Duncan. swift\_readpl is a derivative of SWIFT developed by Jean-Marc Pétit, which swift\_glisse\_ext is based upon.

Usage
-----

### No encounters

In the basic case, integration without close encounters, one must supply an initial particle/planet state and a configuration file. The folder `example/resonance` contains an example state and configuration file which has a large number (200,000) of particles in a semi-major axis range around the 3:1 resonance with Neptune for a duration of 1 million years. The program can be started by running `bin/glisser example/resonance/config.in`. On my computer, this takes approximately 15 minutes to finish. Here are the contents of the configuration file:

    Initial-Time 0
    Time-Step 180
    Final-Time 365e6
    Time-Block-Size 2048
    Track-Interval 2
    Log-Interval 32
    Output-Folder example/resonance-out
    Dump-Interval 20000
    Input-File example/resonance/state.in


*   `Initial-Time`, `Final-Time`, `Time-Step`: exactly what they sound like (see below for comments on units)
*   `Time-Block-Size`: the number of timesteps to integrate at once on the GPU (recommended >500 for GPU-only integrations)
*   `Track-Interval`: the interval (in units of time-blocks) between times where GLISSER should output particle history files
*   `Log-Interval`: the interval (in units of time-blocks) between times where GLISSER should output the status to stdout
*   `Output-Folder`: the folder (absolute or relative) where GLISSER should put all of its output files
*   `Dump-Interval`: the interval (in units of time-blocks) between times where GLISSER should output full binary dumps of particle and planet positions
*   `Input-File`: the input state which contains planet and particle initial positions

Since we are not resolving encounters, particles are removed if they come within 0.5 au of a planet (can be changed by the parameter `Cull-Radius`). The check is made every timestep. In addition, the integration will terminate if any of the planets’ orbits become unbound or fail to converge (e.g. in a situation where two planets come very close to each other).

### Planetary history pre-integration

GLISSER can also be used to only integrate planets to generate a file for planetary history interpolation. `example/threebody-planets` contains an initial condition file with Neptune on a circular orbit and a micro-Jupiter (with Saturn and Uranus removed), which satisfies the restricted three-body problem. The integration lasts for 1000 years, and outputs a planetary history every 100 years. Below is the configuration file.

    Initial-Time 0
    Time-Step 0.2
    Final-Time 1000
    Time-Block-Size 5
    Resolve-Encounters 0
    Track-Interval 1
    Log-Interval 64
    Output-Folder example/threebody-planets-out
    Dump-Interval 0
    Limit-Particle-Count 0
    Input-File example/threebody-planets/state.in
    Swift-History-Interval 100


*   encounter resolution is explicitly disabled by setting `Resolve-Encounters` to 0
*   particle integration on the GPU is disabled by setting `Limit-Particle-Count` to 0
*   planetary history output is enabled by setting `Swift-History-Interval` to 100, which outputs a planetary state for interpolation every 100 timeblocks (100 years)

This takes less than one second on my machine to run, and outputs the file `example/threebody-planets-out/plhist.out` which can be used as an interpolated planetary history.

### Particle integration with close encounters

A planetary history file can be used instead of direct integration of planets. When handling particle close encounters, one must use a planetary history file. `example/threebody-particles` contains a configuration file which uses the planetary history from `example/threebody-planets` to integrate 200,000 particles over 1000 years with close encounters enabled. Below is the configuration file.

    Initial-Time 0
    Time-Step 0.2
    Final-Time 1000
    Time-Block-Size 100
    Resolve-Encounters 1
    Planet-History-File example/threebody-planets-out/plhist.out
    Read-Planet-History 1
    Swift-Path swift/main/swift_glisse_ext
    Track-Interval 1
    Encounter-RH-Factor 3.5
    Cull-Radius 0.5
    Log-Interval 256
    Output-Folder example/threebody-particles-out
    Dump-Interval 20000
    Input-File example/threebody-particles/state.in


*   `Time-Block-Size` should be set such that `Time-Block-Size` \* `Time-Step` is a couple of factors less than the planetary history interval
*   `Resolve-Encounters` is enabled
*   `Planet-History-File` points to the planetary history file
*   `Swift-Path` points to the location of `swift_glisse_ext` on the machine
*   `Track-Interval`, `Log-Interval`, and `Dump-Interval` instead refer to an interval in units of planetary history intervals, instead of timeblocks, when `Resolve-Encounters` is enabled
*   `Encounter-RH-Factor` is the distance in multiples of Hill spheres at which a particle is marked as entering an encounter, and should be set to the same number as `RHSCALE` in `swift/rmvs/rmvs.inc`
*   `Cull-Radius` instead refers to the distance to the sun below which particles should be removed when `Encounter-RH-Factor` is enabled

Inputs
------

### Units

GLISSER uses any natural combination of units where G\=1G = 1G\=1. It is recommended (and in the case of reading from a planetary history file, mandatory) to use the following units:

*   Time: yr
*   Distance: au
*   Mass: 4π2msun4\\pi^2 m\_{sun}4π2msun​ (i.e., the sun has mass 4π24\\pi^24π2)
    Note that the example integration `example/resonance` and the input data file `data/outer-planets.in` use a time unit of days instead, which is not recommended. It is also possible to set `Big-G` in the configuration file to use a different combination of units.

### Input state file

An input state file consists of planetary and particle states at the beginning of the integration. The positions and velocities can be in any inertial frame.

#### ASCII

    Number of planets (including the sun)
    (repeat for each planet, including the sun:)
        mass
        x y z
        vx vy vz
        id
    Number of particles
    (repeat for each particle:)
        x y z
        vx vy vz
        id flags deathtime


#### Binary

    int64   n_planets (including sun)
    repeat n_planets times:
        float64   mass
        int32     id
        float64   x
        float64   y
        float64   z
        float64   vx
        float64   vy
        float64   vz
    int64   n_particles
    repeat n_particles times:
        int32     id
        float64   x
        float64   y
        float64   z
        float64   vx
        float64   vy
        float64   vz
        int16     flags
        float64   deathtime


#### Flags

Particle flags specification:

    High byte  Low byte
    planet id
    <------->
    xxxx xxxx  xxxx xxxx
               ^    ^^^^
              /    / || \
             /    /  | \  Close encounter
           Dead  |   |   Out of bounds
                /    Kepler didn't converge
             Unbound


### Planetary history file

A planetary history file can be used to prove planetary locations instead of a direct integration, if `Read-Planetary-History` is enabled. The planetary history file consists of heliocentric orbital elements, which are interpolated to obtain planetary positions. The format of the planetary history file is as below.

    float64  sun mass
    none192  24 bytes of padding
    repeat until completion:
        float64  time
        int32    n_planets (not including sun)
        none160  20 bytes of padding
        repeat n_planets times:
            int32    planet id
            float32  planet mass
            float32  planet semi-major axis
            float32  planet eccentricity
            float32  planet inclination
            float32  planet longitude of the ascending node
            float32  planet argument of periapsis
            float32  planet mean anomaly


### Configuration

| Name | Description | Default |
| --- | --- | --- |
|`Initial-Time` | The starting time of the integration. | 0 |
| `Time-Step` |The timestep to use in the integration. | 1 |
|`Final-Time`|The time to stop the integration.| |
|`Time-Block-Size`|The timeblock size; the planetary chunk size. The number of timesteps that the GPU will advance in one kernel launch.|1024|
|`Big-G`|The numerical value of G.|1|
|`Encounter-RH-Factor` |**If not 0:**<br> Particles are marked as in an encounter if they come within `Encounter-RH-Factor` \* rHr\_HrH​ of a planet.| 0 |
|`Cull-Radius`|**If `Encounter-RH-Factor` is not set:** <br>Particles are marked as in an encounter if they come within this radius of any planet. <br>**Else:** <br>Particles are marked as in an encounter if they come within this radius of the sun.|0.5|
|`Outer-Limit`|Particles are removed if they are this far away from the sun.|300|
|`Limit-Particle-Count`|Limit the number of particles that are read and integrated from the input file.**If 0:**Performs CPU integration of planets exclusively.|MAX\_INT|
|`Log-Interval` |**If 0:**<br> Progress printing is disabled.<br> **If `Resolve-Encounters` is not set:** <br>Print the current progress every `Log-Interval` timeblocks. <br>**Else:** <br>Print the current progress every `Log-Interval` planetary history intervals.|10 |
|`Track-Interval`|**If 0:**<br> Particle histories are disabled. <br>**If `Resolve-Encounters` is not set:** <br>Write particle history every `Log-Interval` timeblocks.<br> **Else:** <br>Write particle history every `Log-Interval` planetary history intervals.|0|
|`Resync-Interval`|**If `Resolve-Encounters` is not set:** <br> Resort (defragment) the particle arrays every `Resync-Interval` timeblocks.<br> **Else:** <br> This parameter must be 1. |1 |
|`Write-Barycentric-Track` |Write barycentric instead of heliocentric orbital elements to the particle histories if enabled. |1 |
|`Write-RV-Track` |Write position and velocity vectors into track.||
|`Split-Track-File` | **If 0:** <br> Write particle tracks into a single file named `tracks/track.0.out` in the output directory. <br>**Else:** <br>Write particle tracks to files with a maximum size of `Split-Track-File` in bytes, named sequentially in the `tracks` folder in the output directory. |0 |
|`Dump-Interval` |**If 0:** <br> State dumps are disabled. <br>**If `Resolve-Encounters` is not set:** <br> Dump particle and planetary states every `Log-Interval` timeblocks.<br>**Else:**<br> Dump particle nad planetary states every `Log-Interval` planetary history intervals. |1000|
|`Resolve-Encounters` |**If 0:** <br> Particles are removed when coming into an encounter with a planet. <br>**If 1:**<br> Particles are sent to SWIFT when coming into an encounter with a planet. `Resync-Interval` must be 1. `Read-Planetary-History` must be 1. `Swift-Path` must be set. |0|
|`Read-Planetary-History` |Whether to read a planetary history file instead of performing direct planetary integration. |0|
|`Planet-History-File` |The path of the planetary history file. ||
|`Use-Bary-Interpolation` |Whether to use barycentric for orbital interpolation. ||
|`Planet-History-Max-Planet-Count` |The total number of planets that the planetary history file contains.  |16
|`Swift-Path` |The path to swift\_glisse\_ext. ||
|`Swift-Process-Count` |The maximum number of concurrent SWIFT processes to run. |1 |
|`Swift-Process-Min-Particle-Count`|The minimum number of particles on each SWIFT process to warrant launching a new SWIFT process.|10|
|`Swift-Status-Length` |The length of the SWIFT rstat and istat arrays. Must be equal to `NSTAT` in `swift/swift.inc`.  |13
|`Write-Binary-Output` |Whether to write the output state file in binary format. |0 |
|`Read-Binary-Input` |Whether to write the input state file in binary format. |0 |
|`Input-File` |The path of the input state file to read. ||
|`Output-File` |The path of the output folder.||



Output files
------------

### Particle history (track) file

Binary specification:

    repeat:
        float64 time
        int64   n_planets (not including sun)
        repeat n_planets times:
            int32   id
            float32 a
            float32 e
            float32 i
            float32 long. asc. node
            float32 arg. peri
            float32 true anomaly
        int64   n_particles
        repeat n_particles times:
            int32   id
            float32 a
            float32 e
            float32 i
            float32 long. asc. node
            float32 arg. peri
            float32 true anomaly


Particles are ordered by ascending id. Dead particles are not written, so `n_particles` can vary between snapshots. Snapshots occur on regular intervals governed by `Track-Interval`, except if `Read-Planet-History` is set, where outputs can occur on irregular intervals if the planet history intervals are irregular.

Further reading
---------------

Please see also `glisser_manual.md` for a detailed explanation of the code.