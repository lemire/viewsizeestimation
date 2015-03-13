A data warehouse cannot materialize all possible views, hence we must estimate quickly, accurately, and reliably the size of views to determine the best candidates for materialization. Many available techniques for view-size estimation make particular statistical assumptions and their error can be large. Comparatively, unassuming probabilistic techniques are slower, but they estimate accurately and reliability very large view sizes using little memory. We propose five unassuming hashing-based view-size estimation techniques including Stochastic Probabilistic Counting, LogLog Probabilistic Counting, Generalized Counting, Gibbons-Tirthapura, and Adaptive Counting.

More details are available from the paper:

Kamel Aouiche and Daniel Lemire, A Comparison of Five Probabilistic View-Size Estimation Techniques in OLAP, DOLAP 2007, pp. 17-24, 2007.
 http://arxiv.org/abs/cs/0703058.
