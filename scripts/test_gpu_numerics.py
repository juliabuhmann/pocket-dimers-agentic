import jax, jax.numpy as jnp, numpy as np
print(jax.__version__, jax.devices())
rng=np.random.default_rng(0)
def chk(name, f, *a):
    out=np.asarray(jax.jit(f)(*[jnp.asarray(x) for x in a]),dtype=np.float32)
    ref=np.asarray(f(*[np.asarray(x,dtype=np.float32) for x in a]),dtype=np.float32) if False else None
    return out
for dt in [jnp.float32, jnp.bfloat16, jnp.float16]:
    a=rng.standard_normal((8,4,128,64)).astype(np.float32); b=rng.standard_normal((8,4,128,64)).astype(np.float32)
    f=lambda x,y: jnp.einsum('bhqc,bhkc->bhqk',x,y)
    o=np.asarray(jax.jit(f)(jnp.asarray(a,dt),jnp.asarray(b,dt)),dtype=np.float32)
    r=np.einsum('bhqc,bhkc->bhqk',a,b)
    print(jnp.dtype(dt).name,'nan',np.isnan(o).any(),'maxerr',np.abs(o-r).max())
x=jnp.asarray(rng.standard_normal((4,16,300,64)),jnp.float32)
sm=jax.jit(lambda l: jax.nn.softmax(l,axis=-1))(x); print('softmax nan',bool(jnp.isnan(sm).any()))
ln=jax.jit(lambda l:(l-l.mean(-1,keepdims=True))/jnp.sqrt(l.var(-1,keepdims=True)+1e-5))(x); print('ln nan',bool(jnp.isnan(ln).any()))
xs=np.asarray(x)
print('exp ok', np.allclose(np.asarray(jax.jit(jnp.exp)(x)), np.exp(xs), atol=1e-5))
print('sum err', np.abs(np.asarray(jax.jit(lambda l:l.sum(-1))(x))-xs.sum(-1)).max())
print('max err', np.abs(np.asarray(jax.jit(lambda l:l.max(-1))(x))-xs.max(-1)).max())
print('softmax nan frac', float(np.isnan(np.asarray(sm)).mean()))
