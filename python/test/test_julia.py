

def test_julia_bridge(tol=1e-12):
    
    from juliacall import Main as jl
    jl.seval("using Korg")

    v = jl.Korg.air_to_vacuum(4000.0)
    assert abs(v - 4001.1310211436203) < tol
    
