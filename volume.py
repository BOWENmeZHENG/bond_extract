import matplotlib.pyplot as plt
import numpy as np

V = [6.305359999999999445e+02, 7.009642637594585324e+01, 6.907670052949430328e+01, 6.870329976730162969e+01,
     6.844654036320923751e+01, 6.830287951956212567e+01, 6.819326665057194248e+01, 6.807334985883601064e+01,
     6.787037355694791074e+01, 6.780507338582908972e+01, 6.764316500558642531e+01]
N = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]


plt.figure(figsize=(10, 8))
plt.plot(N, np.log(V), '-o', linewidth=2, markersize=10)
plt.xlabel(r"$N \rm _{cycle}$", fontsize=20)
plt.ylabel(r'$L_{\rm box} \rm \, (\AA)$', fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()