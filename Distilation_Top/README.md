## Distilation Top

Distilation top is distilation mixture by reflux. \
We can solve theorical steps using Mccabe -Thiele methode.

In this experiment, we use Ethanol-water solution. And real step is 8.



### Mccabe-Thiele methode(q-line)

1. Plot liquid-gas phase molfraction equilibrium(1 atm) graph. \
   I use NRTL model to find gamma.
2. Plot Rectifying and stripping operating line.
3. Find dot of Rectifying and stripping line same.\
   and plot line between feed's molfractin( $x_D$ )\
   This line is q-line.
   $y=\frac{q}{q-1}x-\frac{x_D}{q-1}$
4. Start at dot of $x_D$-Rectifying line ( $x_D$ is Top solution molfraction)\
   step.

   ![s](image/q-line.png)



### Process Simulation using `DWSIM`

`DWSIM` is open source process simulater.

##### Before reflux

Real experience molefraction data is: 
Feed: 0.2091
Top: 0.6699
Bottom: 0.07046

Simulated data is:

Stage Molar Fractions - Vapor


Stage                              Water             Ethanol


1                                0.42648             0.57352


2                                0.72746             0.27254


3                               0.727454            0.272546


4                               0.727448            0.272552


5                               0.727444            0.272556


6                               0.727439            0.272561


7                               0.727435            0.272565


8                               0.727431            0.272569


9                               0.727427            0.272573


10                              0.949884           0.0501155


Stage 1 and 10 is condnser and reboiler. Net distia

reboiler product molar flwo 0.00407mol/s
