# Generalized-Nonconvex-Approach-for-Low-Tubal-Rank-Tensor-Recovery
The released code of IRTNN algorithm, mainly proposed in the paper "Generalized Nonconvex Approach for Low-Tubal-Rank Tensor Recovery" pubiliahed in TNNLS 2022
-The released IRTNN (also dubbed as IR-t-TNN in the original paper) is the main proposed algorithm in the following paper published in TNNLS-2022 by Hailin Wang et al.:
--URL: https://ieeexplore.ieee.org/abstract/document/9340243
--Citation: 
@ARTICLE{9340243,
  author={Wang, Hailin and Zhang, Feng and Wang, Jianjun and Huang, Tingwen and Huang, Jianwen and Liu, Xinling},
  journal={IEEE Transactions on Neural Networks and Learning Systems}, 
  title={Generalized Nonconvex Approach for Low-Tubal-Rank Tensor Recovery}, 
  year={2022},
  volume={33},
  number={8},
  pages={3305-3319},
  doi={10.1109/TNNLS.2021.3051650}
}

------ Description-----
--1. It is a generalized nonconvex approach for low-tubal-rank tensor recovery (the main paper considers the tensor completion task) based on iterative reweighting strategy

--2. It is suitable for a large family of nonconvex surrogates of the tubal rank of third order tensor under the t-product induced t-SVD framework by Kilmer et al. 

--3. The used tensor singular value nonconvex penalty functions include ETP, Geman, Laplace, Logarithm, Lp, MCP, Capped L1, SCAD and others, as long as the function is continuous, monotonically nondecreasing, and concave function (see Definition 10 in the paper)

--4. The main contributions of this paper are not limited to the algorithm itself, but the provided strict theoretical ananlysis results, where we proved the algorithm converges to a critical point globally with rigorous proofs based on the Kurdyka–Łojasiwicz property, and provide the theoretical guarantees for exact and robust recovery by developing the tensor null space property, showing the advantage of model ability over classic convex method TNN.



--We sincerely recommend the reader the theoretical analysis part, whose proof details are given in the supplement materials: URL: https://ieeexplore.ieee.org/ielx7/5962385/9849214/9340243/supp1-3051650.pdf?arnumber=9340243

--We are very grateful to the authors of the competing methods SNN, TNN, WTNN, PSTNN and IRNN for sharing their codes, which contributes the community a lot.

------May this work is helpful to you!-----

-------------------------------------------wrriten by Hailin Wang, 2022-12-13--------------------------------------
--if you have any question, please contact me without hesitation--email: wanghailin97@163.com
