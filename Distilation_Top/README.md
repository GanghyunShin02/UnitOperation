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

### 📊 Comparison of Distillation Column Simulation and Experimental Results

| Case | Stage (탑 내부 위치) | Water (Molar Fraction) | Ethanol (Molar Fraction) |
| :---: | :---: | :---: | :---: |
| **실제 실험치 (Target)** | **상부 (Top / 1단)** <br> **하부 (Bottom / Bottoms)** | - <br> - | **0.7231** <br> **0.1334** |
| **이론 8단** <br> (단 효율: 0.625) | **1 (Top / Condenser)** <br> 2 <br> 3 <br> 4 <br> 5 <br> 6 <br> 7 <br> 8 <br> 9 <br> **10 (Bottom / Reboiler)** | 0.261422 <br> 0.319773 <br> 0.355627 <br> 0.394237 <br> 0.439235 <br> 0.446717 <br> 0.461929 <br> 0.492531 <br> 0.552562 <br> **0.662485** | **0.738578** <br> 0.680227 <br> 0.644373 <br> 0.605763 <br> 0.560765 <br> 0.553283 <br> 0.538071 <br> 0.507469 <br> 0.447438 <br> **0.337515** |
| **이론 5단** <br> (단 효율: 1.0) | **1 (Top / Condenser)** <br> 2 <br> 3 <br> 4 <br> 5 <br> 6 <br> **7 (Bottom / Reboiler)** | 0.245415 <br> 0.291727 <br> 0.334744 <br> 0.382732 <br> 0.446335 <br> 0.490011 <br> **0.662928** | **0.754585** <br> 0.708273 <br> 0.665256 <br> 0.617268 <br> 0.553665 <br> 0.509989 <br> **0.337072** |
| **이론 7단** <br> (단 효율: 1.0) | **1 (Top / Condenser)** <br> 2 <br> 3 <br> 4 <br> 5 <br> 6 <br> 7 <br> 8 <br> **9 (Bottom / Reboiler)** | 0.228772 <br> 0.264322 <br> 0.294599 <br> 0.323851 <br> 0.355828 <br> 0.395547 <br> 0.449960 <br> 0.492919 <br> **0.662485** | **0.771228** <br> 0.735678 <br> 0.705401 <br> 0.676149 <br> 0.644172 <br> 0.604453 <br> 0.550040 <br> 0.507081 <br> **0.337515** |

We measured top vapor condensered(condensor; stage 1) and bottom(reboiler;stage 10) liquid.
