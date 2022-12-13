# Kerr Comb Simulator

This MATLAB application was developed during my [master thesis](http://urn.kb.se/resolve?urn=urn:nbn:se:liu:diva-188085) project (now also available in docs if link does not work) to aid in conducting numerical simulations of Kerr frequency combs. It provides a graphical user interface to rapidly schedule, execute and inspect simulations of Kerr cavities over a two-dimensional parameter space corresponding to pump frequency detuning and power.

![Example](/docs/Example.PNG)

## Governing equations

The intracavity field $\psi$ can be modelled either by a normalized Lugiato-Lefever Equation (LLE) 

$$ \frac{\partial \psi}{\partial \eta} = -(\alpha + i \delta) \psi - i s \frac{\partial^2 \psi}{\partial \tau^2} + i |\psi|^2 \psi + \sqrt{\alpha} \psi_{\text{in}},$$

where $\eta$ is the slow time (evolution variable) and $\tau$ the fast time (field profile), or by a normalized Ikeda map

$$ \frac{\partial \psi^m}{\partial \xi} = -\frac{\alpha}{2} \psi^m - i s \frac{\partial^2 \psi^m}{\partial \tau^2} + i |\psi^m|^2 \psi^m,$$

$$\psi^{m + 1}(0, \tau) = \sqrt{1 - \alpha} e^{-i \delta} \psi^m(1, \tau) + \sqrt{\alpha} \psi_{\text{in}}, $$

where the superscript of the field indicates the current roundtrip and the second equation constitute the time-varying boundary condition coupling the fields from successive roundtrips. The equations are expressed using the normalized variables

$$ \xi = \frac{z}{L}, \quad \tau = \sqrt{\frac{2}{|\beta_2| L}} \ t, \quad \psi = \sqrt{L \gamma} \ A, $$

where $z$ is the longitudinal distance in the cavity, $L$ is the cavity length, $t$ is the fast time, $\beta_2 = s |\beta_2|$ is the group velocity dispersion at the pump frequency $\omega_0$, $A$ is the pulse envelope amplitude (scaled such that $|A|^2$ is measured in units of power) and $\gamma = \omega_0 n_2 / c A_{\text{eff}}$ is the nonlinear coefficient with $n_2$ being the Kerr coefficient, $c$ the vacuum speed of light and $A_{\text{eff}}$ the effective mode area. Further, $\alpha$ is the total cavity loss related to the cavity finesse through $\mathcal{F} = \pi / \alpha$ which equals the power coupling coefficient since critical coupling is assumed, $\delta$ is the pump frequency detuning and $\psi_{\text{in}} = \sqrt{L \gamma} \ A_{\text{in}}$ is the normalized pump field.

The two models were normalized with the intent to reduce the parameter space to two dimensions (after choosing the resonator finesse and group velocity dispersion sign), where both models use the same parameter space. This was to make the produced phase diagrams directly comparable since the project was in part intended to look at the differences between the often-employed LLE model and the more general Ikeda map model.

## Basic usage guide

A batch, representing a collection of systems, can be created or loaded from the 'File' tab. All systems in a batch share certain parameters (resonator finesse and group-velocity dispersion sign) and the model used to simulate them (LLE or Ikeda map), along with a few other tuning parameters (see tooltips for documentation).

After a batch is created it can be populated using the 'Create system' panel. Some of the options include initial field (noise, Gaussian, sech), the size of the circular buffer storing the field history for playback, the initial and final values for the pump detuning and power (can be set graphically by enabling the 'Fig' toggles).

When a system is created it is automatically selected, highlighting it on the parameter map and in the system table. When one or more system is selected, the 'Edit system' and 'Fork system' panels replace the 'Display' and 'Create system' ones. From here, collective actions on selected systems can be performed, such as renewing budgets and 'forking' systems using the fields of selected systems as initial values. To select other systems, either use the 'Select' button and brush systems in the parameter map or (shift) click systems in the table. To clear the selection, select click the map or table.

Once one or more systems have been scheduled for simulation, press the 'Run' button to execute them. To execute simulations in parallel, check the 'Parallel execution' box (requires Parallel computing toolbox). When a system is finished executing it reports the expended budget to the progress bar.

When a system finish simulating, it becomes available for inspection when selected. The right-side windows enable playback of the field, spectrum and CW reference evolution through field snapshots captured each series (as far back as the buffer spans). Additionally, the convergence of the system (field variation between series) can be viewed over the entire field history with accompanying color indicating the automatic classification (sweep, stationary, fluctuating, CW). Check the 'Legend' box for color correspondence.


## Compatibility

The app has been tested on MATLAB R2021b using Windows 10. There are some bugs that I know about and likely far more that I don't. 

In order to enable graphical selection ('brushing') on the parameter map, the use of the undocumented property 'BrushData' was necessary. In the (unlikely) scenario that this property is changed in future MATLAB releases, this might break things.
