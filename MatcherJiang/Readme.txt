MatcherJiang.

=====================================================================================================

Jiang and Yau (2000). Jiang X. and Yau W.Y., “Fingerprint Minutiae Matching Based on the Local and
Global Structures,” in Proc. Int. Conf. on Pattern Recognition (15th), vol. 2, pp. 1042–1045, 2000.
=====================================================================================================

A local fingerprint matching algorithm. It works by computing the nearest neighbours minutiae for each minutia
and process them in three phases:

Phase I) Feature Vector Computation: A set of feature vector for each minutiae is computed. It is formed bu distance,
		angles, directins, ridge counts and types of minutiae.

Phase II) Strict Matching: Solely minutiae are paired. The best set of matching minutiae are selected.

Phase III) Consolidation: A single transformation is done for every minutiae considering the best minutiae (template
		and input) as polar coordinated. Tolerance boxes are considered and the matching score is computed.
