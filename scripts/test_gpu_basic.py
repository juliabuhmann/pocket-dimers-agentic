import jax, jax.numpy as jnp, numpy as np
rng=np.random.default_rng(0); x=rng.standard_normal((4,16,300,64)).astype(np.float32)
d=jnp.asarray(x); print('roundtrip', np.array_equal(np.asarray(d),x))
for n,f,r in [('x+1',lambda a:a+1,x+1),('x*x',lambda a:a*a,x*x),('exp',jnp.exp,np.exp(x)),('sum',lambda a:a.sum(-1),x.sum(-1)),('tanh',jnp.tanh,np.tanh(x))]:
    o=np.asarray(jax.jit(f)(d)); print(n,'maxerr',np.abs(o-r).max(),'nan',np.isnan(o).any())
