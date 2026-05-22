# Distilation Top

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


5. Calculated stage efficent using computer is 0.625 but hand used is 0.875.


### Process Simulation using `DWSIM`

`DWSIM` is open source process simulater.

##### Before reflux

Real experience molefraction data is: 
Feed: 0.2091
Top: 0.6699
Bottom: 0.07046

Simulated data 

Reflux ratio 0.2 (Realy in top, should be reflux)\
Except feed stage, condensor, reboiler, stage efficent is 0.625.

##### Stage Conditions & Flows

| Stage | P (Pa) | T (K) | V (mol/s) | L (mol/s) | LSS (mol/s) |
|---|---:|---:|---:|---:|---:|
| 1 | 101325 | 352.358 | 0 | 0.00108394 | 0.00541969 |
| 2 | 101325 | 353.721 | 0.00650363 | 0.00105401 | 0 |
| 3 | 101325 | 354.863 | 0.00647371 | 0.00103389 | 0 |
| 4 | 101325 | 355.642 | 0.00645359 | 0.00102349 | 0 |
| 5 | 101325 | 356.089 | 0.00644318 | 0.0184003 | 0 |
| 6 | 101326 | 356.306 | 0.00665356 | 0.0183749 | 0 |
| 7 | 101326 | 356.763 | 0.00662817 | 0.0183277 | 0 |
| 8 | 101326 | 357.771 | 0.00658098 | 0.0182485 | 0 |
| 9 | 101326 | 360.043 | 0.00650172 | 0.0181347 | 0 |
| 10 | 101326 | 367.005 | 0.00638797 | 0.0117467 | 0 |

##### Stage Molar Fractions - Vapor

| Stage | Water | Ethanol |
|---|---:|---:|
| 1 | 0.302181 | 0.697819 |
| 2 | 0.400089 | 0.599911 |
| 3 | 0.432703 | 0.567297 |
| 4 | 0.450720 | 0.549280 |
| 5 | 0.458863 | 0.541137 |
| 6 | 0.481298 | 0.518702 |
| 7 | 0.508354 | 0.491646 |
| 8 | 0.559018 | 0.440982 |
| 9 | 0.647599 | 0.352401 |
| 10 | 0.779132 | 0.220868 |

##### Stage Molar Fractions - Liquid 1

| Stage | Water | Ethanol |
|---|---:|---:|
| 1 | 0.400089 | 0.599911 |
| 2 | 0.600403 | 0.399597 |
| 3 | 0.716129 | 0.283871 |
| 4 | 0.770084 | 0.229916 |
| 5 | 0.794058 | 0.205942 |
| 6 | 0.804250 | 0.195750 |
| 7 | 0.823204 | 0.176796 |
| 8 | 0.855912 | 0.144088 |
| 9 | 0.903551 | 0.0964486 |
| 10 | 0.971212 | 0.0287885 |


We measured top vapor condensered(condensor; stage 1) and bottom(reboiler;stage 10) liquid.
