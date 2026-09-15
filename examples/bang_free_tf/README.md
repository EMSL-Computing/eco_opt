# Bang-bang minimum-time double integrator

Free-final-time rest-to-rest transfer with asymmetric control bounds.
The model returns a bare `dx` (no signal tree).

```
minimize   tf
subject to xdot = v, vdot = u,  -2 <= u <= 1
           (x,v)(0) = (0,0),  (x,v)(tf) = (300,0)
```

Analytic bang-bang: switch from `u=+1` to `u=-2` at `t=20`, `tf*=30`.

## Running

```bash
cd examples/bang_free_tf
python bang_free_tf.py
python plot_bang_free_tf.py
python -m pytest ../../tests/test_bang_free_tf.py -q
```
