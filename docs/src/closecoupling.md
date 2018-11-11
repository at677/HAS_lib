
# Close coupling formalism

Upercase vectors denote vectors parralel to the surface. For example ``\vec R,
\vec K_i``

Subsitute the potential \$V\$ and the wavefunction \$\psi\$ as ourier series time-independend Schrödinger equation. The potential is written as:

```math
V(\vec r) = \sum_G V_G(z) e^{i \vec G \cdot \vec R}
```

and the wavefunction can be written as:

```math
\psi(\vec r) = \sum_G \psi_G(z) e^{i (\vec G + \vec K_i) \cdot \vec R} 
```

The timeindependend Schrödinger equation for the scattering Problem has the
following form.

```math
\left[-\Delta - \vec k_i^2 + V(\vec r)\right] \psi(\vec r) = 0
```

The first part ``(-\Delta \psi(\vec r))`` can be rewritten using the partial
derivatives and the fourier series of ``\psi``.

```math
\begin{align}
&-(\partial^2_x + \partial^2_y + \partial^2_z) \sum_G \psi_G(z) e^{i (\vec G + \vec K_i) \cdot \vec R} =\\

&- \sum_G e^{i (\vec G + \vec K_i) \cdot \vec R} \partial^2_z \psi_G(z)  -
\sum_G \psi_G(z) (\partial^2_x+\partial^2_y) e^{i (\vec G + \vec K_i) \cdot \vec R} =\\

&- \sum_G e^{i (\vec G + \vec K_i) \cdot \vec R} \psi_G^{''}(z)  +
\sum_G \psi_G(z) \left((G_x + K_{i_x})^2+(G_y + K_{i_y})^2\right) e^{i (\vec G + \vec K_i) \cdot \vec R}=\\

&- \sum_G e^{i (\vec G + \vec K_i) \cdot \vec R} \psi_G^{''}(z)  +
\sum_G \psi_G(z) (\vec G + \vec K_i)^2 e^{i (\vec G + \vec K_i) \cdot \vec
R}=\\

&\sum_G e^{i (\vec G + \vec K_i) \cdot \vec R} \left((\vec G + \vec K_i)^2 \psi_G(z) - \psi_G^{''}(z)\right)

\end{align}
```

Together with the second part ``(\vec k_i^2  \psi(\vec r))`` the equations looks
like this:

```math
\sum_G e^{i (\vec G + \vec K_i) \cdot \vec R} \left[\left((\vec G + \vec K_i)^2-\vec k_i^2\right) \psi_G(z) - \psi_G^{''}(z)\right]
```

The third part (``V(\vec r) \psi(\vec r)``) can be written as:

```math
\begin{align}
&\sum_G V_G(z) e^{i \vec G \cdot \vec R} \sum_G \psi_G(z) e^{i (\vec G + \vec K_i) \cdot \vec R} = \\
&\sum_{\vec G,\vec G'} V_G(z) e^{i (\vec G + \vec G' + \vec K_i) \cdot \vec R}\psi_{G'}(z) 
\end{align}
```

Combining all three parts and multiplying with ``e^{-i (\vec G''+\vec K_i) \cdot \vec R}``
results in:
```math
\sum_G e^{i (\vec G - \vec G^{''}) \cdot \vec R} \left[\left((\vec G + \vec K_i)^2-\vec k_i^2\right) \psi_G(z) - \psi_G^{''}(z)\right]+\\
\sum_{\vec G,\vec G'} V_G(z) e^{i (\vec G + \vec G' - \vec G^{''}) \cdot \vec R}\psi_{G'}(z) = 0
```

Integrating over der unit cell will eliminate the exponential terms by
introducing dirac deltas.

```math
\sum_G \int_U e^{i (\vec G - \vec G^{''}) \cdot \vec R} d\vec R \left[\left((\vec G + \vec K_i)^2-\vec k_i^2\right) \psi_G(z) - \psi_G^{''}(z)\right]+\\
\sum_{\vec G,\vec G'} V_G(z)  \psi_{G'}(z)\int_U e^{i (\vec G + \vec G' - \vec G^{''}) \cdot \vec R} d\vec R = \\

\sum_G A_U\delta(\vec G - \vec G^{''}) \left[\left((\vec G + \vec K_i)^2-\vec k_i^2\right) \psi_G(z) - \psi_G^{''}(z)\right]+\\
\sum_{\vec G,\vec G'} V_G(z)  \psi_{G'}(z)A_U\delta(\vec G + \vec G' - \vec G^{''}) = \\

A_U \left[\left((\vec G^{''} + \vec K_i)^2-\vec k_i^2\right) \psi_{G^{''}}(z) - \psi_{G^{''}}^{''}(z)\right]+
\sum_{\vec G'} V_{G^{''}-G^{'}}(z)  \psi_{G^{'}}(z)A_U 
```

With this the Schrödinger equation transforms to:
```math
\psi_{G^{''}}^{''}(z) 
 = \sum_{\vec G'} V_{G^{''}-G^{'}}(z)  \psi_{G^{'}}(z)
- \left(\vec k_i^2 - (\vec G^{''} + \vec K_i)^2\right) \psi_{G^{''}}(z) 
```


