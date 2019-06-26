---


---

<h1 id="glisser---gpu-long-term-integrator-for-solar-system-evolution-regularized">GLISSER - GPU Long-term Integrator for Solar System Evolution (Regularized)</h1>
<h2 id="introduction">Introduction</h2>
<p>GLISSER is an integrator which makes use of a GPU card to perform long-term and long-scale integrations of massless test particles in a solar system.</p>
<h3 id="functions">Functions</h3>
<p>GLISSER is capable of:</p>
<ul>
<li>GPU-only integration using either an interpolated planetary history file or planetary initial conditions, which does not support particle close encounters</li>
<li>GPU-CPU hybrid integration using an interpolated planetary history file, which supports particle close encounters</li>
</ul>
<h3 id="requirements">Requirements</h3>
<ul>
<li>A CUDA-enabled GPU (CUDA 8 recommended)</li>
<li>Linux</li>
</ul>
<h2 id="installation">Installation</h2>
<h3 id="compiling-the-main-code">Compiling the main code</h3>
<p>GLISSER functions as a standalone code. Simply type <code>make</code> in the project root, which will generate the executable <code>bin/glisse</code>. It may be necessary to edit the variable <code>CUDA_ARCH</code> in <code>Makefile</code> to match your system’s GPU architecture.</p>
<h3 id="compiling-the-addon-swift-code">Compiling the addon SWIFT code</h3>
<p>Navigate into <code>swift/</code> and run <code>./CompileAll.sh</code>, which generates the executable <code>swift/main/swift_glisse_ext</code>. It may be necessary to edit the variables <code>NPLMAX</code> and <code>NTPMAX</code> in <code>swift/main/swift.inc</code> to increase the maximum number of particles or planets that SWIFT can integrate at once.</p>
<blockquote>
<p><a href="https://www.boulder.swri.edu/~hal/swift.html">SWIFT</a> is developed by Harold Levison  and Martin Duncan. swift_readpl is a derivative of SWIFT developed by Jean-Marc Pétit, which swift_glisse_ext is based upon.</p>
</blockquote>
<h2 id="usage">Usage</h2>
<h3 id="no-encounters">No encounters</h3>
<p>In the basic case, integration without close encounters, one must supply an initial particle/planet state and a configuration file. The folder <code>example/resonance</code> contains an example state and configuration file which has a large number (200,000) of particles in a semi-major axis range around the 3:1 resonance with Neptune for a duration of 1 million years. The program can be started by running <code>bin/glisser example/resonance/config.in</code>. On my computer, this takes approximately 15 minutes to finish. Here are the contents of the configuration file:</p>
<pre><code>Initial-Time 0
Time-Step 180
Final-Time 365e6
Time-Block-Size 2048
Track-Interval 2
Log-Interval 32
Output-Folder example/resonance-out
Dump-Interval 20000
Input-File example/resonance/state.in
</code></pre>
<ul>
<li><code>Initial-Time</code>, <code>Final-Time</code>, <code>Time-Step</code>: exactly what they sound like (see below for comments on units)</li>
<li><code>Time-Block-Size</code>: the number of timesteps to integrate at once on the GPU (recommended &gt;500 for GPU-only integrations)</li>
<li><code>Track-Interval</code>: the interval (in units of time-blocks) between times where GLISSER should output particle history files</li>
<li><code>Log-Interval</code>: the interval (in units of time-blocks) between times where GLISSER should output the status to stdout</li>
<li><code>Output-Folder</code>: the folder (absolute or relative) where GLISSER should put all of its output files</li>
<li><code>Dump-Interval</code>: the interval (in units of time-blocks) between times where GLISSER should output full binary dumps of particle and planet positions</li>
<li><code>Input-File</code>: the input state which contains planet and particle initial positions</li>
</ul>
<p>Since we are not resolving encounters, particles are removed if they come within 0.5 au of a planet (can be changed by the parameter <code>Cull-Radius</code>). The check is made every timestep. In addition, the integration will terminate if any of the planets’ orbits become unbound or fail to converge (e.g. in a situation where two planets come very close to each other).</p>
<h3 id="planetary-history-pre-integration">Planetary history pre-integration</h3>
<p>GLISSER can also be used to only integrate planets to generate a file for planetary history interpolation. <code>example/threebody-planets</code> contains an initial condition file with Neptune on a circular orbit and a micro-Jupiter (with Saturn and Uranus removed), which satisfies the restricted three-body problem. The integration lasts for 1000 years, and outputs a planetary history every 100 years. Below is the configuration file.</p>
<pre><code>Initial-Time 0
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
</code></pre>
<ul>
<li>encounter resolution is explicitly disabled by setting <code>Resolve-Encounters</code> to 0</li>
<li>particle integration on the GPU is disabled by setting <code>Limit-Particle-Count</code> to 0</li>
<li>planetary history output is enabled by setting <code>Swift-History-Interval</code> to 100, which outputs a planetary state for interpolation every 100 timeblocks (100 years)</li>
</ul>
<p>This takes less than one second on my machine to run, and outputs the file <code>example/threebody-planets-out/plhist.out</code> which can be used as an interpolated planetary history.</p>
<h3 id="particle-integration-with-close-encounters">Particle integration with close encounters</h3>
<p>A planetary history file can be used instead of direct integration of planets. When handling particle close encounters, one must use a planetary history file. <code>example/threebody-particles</code> contains a configuration file which uses the planetary history from <code>example/threebody-planets</code> to integrate 200,000 particles over 1000 years with close encounters enabled. Below is the configuration file.</p>
<pre><code>Initial-Time 0
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
</code></pre>
<ul>
<li><code>Time-Block-Size</code> should be set such that <code>Time-Block-Size</code> * <code>Time-Step</code> is a couple of factors less than the planetary history interval</li>
<li><code>Resolve-Encounters</code> is enabled</li>
<li><code>Planet-History-File</code> points to the planetary history file</li>
<li><code>Swift-Path</code> points to the location of <code>swift_glisse_ext</code> on the machine</li>
<li><code>Track-Interval</code>, <code>Log-Interval</code>, and <code>Dump-Interval</code> instead refer to an interval in units of planetary history intervals, instead of timeblocks, when <code>Resolve-Encounters</code> is enabled</li>
<li><code>Encounter-RH-Factor</code> is the distance in multiples of Hill spheres at which a particle is marked as entering an encounter, and should be set to the same number as <code>RHSCALE</code> in <code>swift/rmvs/rmvs.inc</code></li>
<li><code>Cull-Radius</code> instead refers to the distance to the sun below which particles should be removed when <code>Encounter-RH-Factor</code> is enabled</li>
</ul>
<h2 id="inputs">Inputs</h2>
<h3 id="units">Units</h3>
<p>GLISSER uses any natural combination of units where <span class="katex--inline"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><mi>G</mi><mo>=</mo><mn>1</mn></mrow><annotation encoding="application/x-tex">G = 1</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="base"><span class="strut" style="height: 0.68333em; vertical-align: 0em;"></span><span class="mord mathit">G</span><span class="mspace" style="margin-right: 0.277778em;"></span><span class="mrel">=</span><span class="mspace" style="margin-right: 0.277778em;"></span></span><span class="base"><span class="strut" style="height: 0.64444em; vertical-align: 0em;"></span><span class="mord">1</span></span></span></span></span>. It is recommended (and in the case of reading from a planetary history file, mandatory) to use the following units:</p>
<ul>
<li>Time: yr</li>
<li>Distance: au</li>
<li>Mass: <span class="katex--inline"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><mn>4</mn><msup><mi>π</mi><mn>2</mn></msup><msub><mi>m</mi><mrow><mi>s</mi><mi>u</mi><mi>n</mi></mrow></msub></mrow><annotation encoding="application/x-tex">4\pi^2 m_{sun}</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="base"><span class="strut" style="height: 0.964108em; vertical-align: -0.15em;"></span><span class="mord">4</span><span class="mord"><span class="mord mathit" style="margin-right: 0.03588em;">π</span><span class="msupsub"><span class="vlist-t"><span class="vlist-r"><span class="vlist" style="height: 0.814108em;"><span class="" style="top: -3.063em; margin-right: 0.05em;"><span class="pstrut" style="height: 2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight">2</span></span></span></span></span></span></span></span><span class="mord"><span class="mord mathit">m</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height: 0.151392em;"><span class="" style="top: -2.55em; margin-left: 0em; margin-right: 0.05em;"><span class="pstrut" style="height: 2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mathit mtight">s</span><span class="mord mathit mtight">u</span><span class="mord mathit mtight">n</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height: 0.15em;"><span class=""></span></span></span></span></span></span></span></span></span></span> (i.e., the sun has mass <span class="katex--inline"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><mn>4</mn><msup><mi>π</mi><mn>2</mn></msup></mrow><annotation encoding="application/x-tex">4\pi^2</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="base"><span class="strut" style="height: 0.814108em; vertical-align: 0em;"></span><span class="mord">4</span><span class="mord"><span class="mord mathit" style="margin-right: 0.03588em;">π</span><span class="msupsub"><span class="vlist-t"><span class="vlist-r"><span class="vlist" style="height: 0.814108em;"><span class="" style="top: -3.063em; margin-right: 0.05em;"><span class="pstrut" style="height: 2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight">2</span></span></span></span></span></span></span></span></span></span></span></span>)<br>
Note that the example integration <code>example/resonance</code> and the input data file <code>data/outer-planets.in</code> use a time unit of days instead, which is not recommended. It is also possible to set <code>Big-G</code> in the configuration file to use a different combination of units.</li>
</ul>
<h3 id="input-state-file">Input state file</h3>
<p>An input state file consists of planetary and particle states at the beginning of the integration. The positions and velocities can be in any inertial frame.</p>
<h4 id="ascii">ASCII</h4>
<pre><code>Number of planets (including the sun)
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
</code></pre>
<h4 id="binary">Binary</h4>
<pre><code>int64   n_planets (including sun)
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
</code></pre>
<h4 id="flags">Flags</h4>
<p>Particle flags specification:</p>
<pre><code>High byte  Low byte
planet id
&lt;-------&gt; 
xxxx xxxx  xxxx xxxx
           ^    ^^^^
          /    / || \
         /    /  | \  Close encounter   
       Dead  |   |   Out of bounds
            /    Kepler didn't converge
         Unbound
</code></pre>
<h3 id="planetary-history-file">Planetary history file</h3>
<p>A planetary history file can be used to prove planetary locations instead of a direct integration, if <code>Read-Planetary-History</code> is enabled. The planetary history file consists of heliocentric orbital elements, which are interpolated to obtain planetary positions. The format of the planetary history file is as below.</p>
<pre><code>float64  sun mass
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
</code></pre>
<h3 id="configuration">Configuration</h3>

<table>
<thead>
<tr>
<th>Name</th>
<th>Description</th>
<th>Default</th>
</tr>
</thead>
<tbody>
<tr>
<td><code>Initial-Time</code></td>
<td>The starting time of the integration.</td>
<td>0</td>
</tr>
<tr>
<td><code>Time-Step</code></td>
<td>The timestep to use in the integration.</td>
<td>1</td>
</tr>
<tr>
<td><code>Final-Time</code></td>
<td>The time to stop the integration.</td>
<td></td>
</tr>
<tr>
<td><code>Time-Block-Size</code></td>
<td>The timeblock size; the planetary chunk size. The number of timesteps that the GPU will advance in one kernel launch.</td>
<td>1024</td>
</tr>
<tr>
<td><code>Big-G</code></td>
<td>The numerical value of G.</td>
<td>1</td>
</tr>
<tr>
<td><code>Encounter-RH-Factor</code></td>
<td><strong>If not 0:</strong> <br> Particles are marked as in an encounter if they come within <code>Encounter-RH-Factor</code> * <span class="katex--inline"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><msub><mi>r</mi><mi>H</mi></msub></mrow><annotation encoding="application/x-tex">r_H</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="base"><span class="strut" style="height: 0.58056em; vertical-align: -0.15em;"></span><span class="mord"><span class="mord mathit" style="margin-right: 0.02778em;">r</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height: 0.328331em;"><span class="" style="top: -2.55em; margin-left: -0.02778em; margin-right: 0.05em;"><span class="pstrut" style="height: 2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mathit mtight" style="margin-right: 0.08125em;">H</span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height: 0.15em;"><span class=""></span></span></span></span></span></span></span></span></span></span> of a planet.</td>
<td>0</td>
</tr>
<tr>
<td><code>Cull-Radius</code></td>
<td><strong>If <code>Encounter-RH-Factor</code> is not set:</strong> <br> Particles are marked as in an encounter if they come within this radius of any planet. <br> <strong>Else:</strong> <br> Particles are marked as in an encounter if they come within this radius of the sun.</td>
<td>0.5</td>
</tr>
<tr>
<td><code>Outer-Limit</code></td>
<td>Particles are removed if they are this far away from the sun.</td>
<td>300</td>
</tr>
<tr>
<td><code>Limit-Particle-Count</code></td>
<td>Limit the number of particles that are read and integrated from the input file. <br> <strong>If 0:</strong> <br> Performs CPU integration of planets exclusively.</td>
<td>MAX_INT</td>
</tr>
<tr>
<td><code>Log-Interval</code></td>
<td><strong>If 0:</strong> <br> Progress printing is disabled. <br> <strong>If <code>Resolve-Encounters</code> is not set:</strong> <br> Print the current progress every <code>Log-Interval</code> timeblocks. <br> <strong>Else:</strong> <br> Print the current progress every <code>Log-Interval</code> planetary history intervals.</td>
<td>10</td>
</tr>
<tr>
<td><code>Track-Interval</code></td>
<td><strong>If 0:</strong> <br> Particle histories are disabled. <br> <strong>If <code>Resolve-Encounters</code> is not set:</strong> <br> Write particle history every <code>Log-Interval</code> timeblocks. <br> <strong>Else:</strong> <br> Write particle history every <code>Log-Interval</code> planetary history intervals.</td>
<td>0</td>
</tr>
<tr>
<td><code>Resync-Interval</code></td>
<td><strong>If <code>Resolve-Encounters</code> is not set:</strong> <br> Resort (defragment) the particle arrays every <code>Resync-Interval</code> timeblocks. <br> <strong>Else:</strong> <br> This parameter must be 1.</td>
<td>1</td>
</tr>
<tr>
<td><code>Write-Barycentric-Track</code></td>
<td>Write barycentric instead of heliocentric orbital elements to the particle histories if enabled.</td>
<td>1</td>
</tr>
<tr>
<td><code>Split-Track-File</code></td>
<td><strong>If 0:</strong> <br> Write particle tracks into a single file named <code>tracks/track.0.out</code> in the output directory. <br> <strong>Else:</strong> Write particle tracks to files with a maximum size of <code>Split-Track-File</code> in bytes, named sequentially in the <code>tracks</code> folder in the output directory.</td>
<td>0</td>
</tr>
<tr>
<td><code>Dump-Interval</code></td>
<td><strong>If 0:</strong> <br> State dumps are disabled. <br> <strong>If <code>Resolve-Encounters</code> is not set:</strong> <br> Dump particle and planetary states every <code>Log-Interval</code> timeblocks. <br> <strong>Else:</strong> <br> Dump particle nad planetary states every <code>Log-Interval</code> planetary history intervals.</td>
<td>1000</td>
</tr>
<tr>
<td><code>Resolve-Encounters</code></td>
<td><strong>If 0:</strong> <br> Particles are removed when coming into an encounter with a planet. <br> <strong>If 1:</strong> <br> Particles are sent to SWIFT when coming into an encounter with a planet. <code>Resync-Interval</code> must be 1. <code>Read-Planetary-History</code> must be 1. <code>Swift-Path</code> must be set.</td>
<td>0</td>
</tr>
<tr>
<td><code>Read-Planetary-History</code></td>
<td>Whether to read a planetary history file instead of performing direct planetary integration.</td>
<td>0</td>
</tr>
<tr>
<td><code>Planet-History-File</code></td>
<td>The path of the planetary history file.</td>
<td></td>
</tr>
<tr>
<td><code>Planet-History-Max-Planet-Count</code></td>
<td>The total number of planets that the planetary history file contains.</td>
<td>16</td>
</tr>
<tr>
<td><code>Swift-Path</code></td>
<td>The path to swift_glisse_ext.</td>
<td></td>
</tr>
<tr>
<td><code>Swift-Process-Count</code></td>
<td>The maximum number of concurrent SWIFT processes to run.</td>
<td>1</td>
</tr>
<tr>
<td><code>Swift-Process-Min-Particle-Count</code></td>
<td>The minimum number of particles on each SWIFT process to warrant launching a new SWIFT process.</td>
<td>10</td>
</tr>
<tr>
<td><code>Swift-Status-Length</code></td>
<td>The length of the SWIFT rstat and istat arrays. Must be equal to <code>NSTAT</code> in <code>swift/swift.inc</code>.</td>
<td>13</td>
</tr>
<tr>
<td><code>Write-Binary-Output</code></td>
<td>Whether to write the output state file in binary format.</td>
<td>0</td>
</tr>
<tr>
<td><code>Read-Binary-Input</code></td>
<td>Whether to write the input state file in binary format.</td>
<td>0</td>
</tr>
<tr>
<td><code>Input-File</code></td>
<td>The path of the input state file to read.</td>
<td></td>
</tr>
<tr>
<td><code>Output-File</code></td>
<td>The path of the output folder.</td>
<td></td>
</tr>
</tbody>
</table><h2 id="output-files">Output files</h2>
<h3 id="particle-history-track-file">Particle history (track) file</h3>
<p>Binary specification:</p>
<pre><code>repeat:
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
</code></pre>
<p>Particles are ordered by ascending id. Dead particles are not written, so <code>n_particles</code> can vary between snapshots. Snapshots occur on regular intervals governed by <code>Track-Interval</code>, except if <code>Read-Planet-History</code> is set, where outputs can occur on irregular intervals if the planet history intervals are irregular.</p>
<h2 id="further-reading">Further reading</h2>
<p>Please see also <code>glisser_manual.md</code> for a detailed explanation of the code.</p>

