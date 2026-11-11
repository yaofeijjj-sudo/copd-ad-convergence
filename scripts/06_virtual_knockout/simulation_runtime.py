"""Keep the out-of-distribution table that CellOracle's simulate_shift computes internally.

The wrapper calls CellOracle's own method and returns its result unchanged, so the table can be
exported without evaluating it a second time.
"""


def simulate_with_ood_diagnostic(model, **kwargs):
    name = "evaluate_simulated_gene_distribution_range"
    original = getattr(model, name)
    had_instance_attribute = name in model.__dict__
    retained = []

    def retain_result(*args, **kw):
        result = original(*args, **kw)
        retained.append(result)
        return result

    setattr(model, name, retain_result)
    try:
        model.simulate_shift(**kwargs)
    finally:
        if had_instance_attribute:
            setattr(model, name, original)
        else:
            delattr(model, name)
    assert len(retained) == 1, "simulate_shift must evaluate OOD once"
    return retained[0].copy(deep=True)
