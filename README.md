# Hybrid-Field Channel Estimation for XL-MIMO Systems with Stochastic Gradient Pursuit Algorithm
This is the code for [1]. 

[1] Hao Lei, Jiayi Zhang, Zhe Wang, Bo Ai, and Derrick Wing Kwan Ng, “Hybrid-Field Channel Estimation for XL-MIMO Systems with Stochastic Gradient Pursuit Algorithm”, 
  IEEE Transactions on Signal Processing, to appear, 2024.


# Abstract
Extremely large-scale multiple-input multipleoutput (XL-MIMO) is crucial for satisfying the high data rate requirements of the sixth-generation (6G) wireless networks.
In this context, ensuring accurate acquisition of channel state information (CSI) with low complexity becomes imperative. Moreover, deploying an extremely large antenna array at the
base station (BS) might result in some scatterers being located in near-field, while others are situated in far-field, leading
to a hybrid-field communication scenario. To address these challenges, this paper introduces two stochastic gradient pursuit
(SGP)-based schemes for the hybrid-field channel estimation in two scenarios. For the first scenario in which the prior
knowledge of the specific proportion of the number of near-field and far-field channel paths is known, the scheme can effectively leverage the angular-domain sparsity of the far-field channels
and the polar-domain sparsity of the near-field channels such that the channel estimation in these two fields can be performed separately. For the second scenario which the proportion
is not available, we propose an off-grid SGP-based channel estimation scheme, which iterates through the values of the proportion parameter based on a criterion before performing
the hybrid-field channel estimation. We demonstrate numerically that both of the proposed channel estimation schemes achieve superior performance in terms of both estimation accuracy
and achievable rates while enjoying lower computational complexity compared with existing schemes. Additionally, we reveal that as the number of antennas at the UE increases, the
normalized mean square error (NMSE) performances of the proposed schemes remain basically unchanged, while the NMSE performances of existing ones improve. Remarkably, even in
this scenario, the proposed schemes continue to outperform the existing ones.





# License and Referencing
If you in any way use this code for research that results in publications, please cite our original article listed above ([1]).
[1] Hao Lei, Jiayi Zhang, Zhe Wang, Bo Ai, and Derrick Wing Kwan Ng, “Hybrid-Field Channel Estimation for XL-MIMO Systems with Stochastic Gradient Pursuit Algorithm”, 
  IEEE Transactions on Signal Processing, to appear, 2024.
