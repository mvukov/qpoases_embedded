""" Defines rule for a test data generator.
"""

def _data_generator_impl(ctx):
    output = ctx.actions.declare_file(ctx.attr.output)
    args = ctx.actions.args()
    args.add(output)
    ctx.actions.run(
        outputs = [output],
        executable = ctx.executable.generator,
        arguments = [args],
    )
    return [DefaultInfo(files = depset([output]))]

data_generator = rule(
    implementation = _data_generator_impl,
    attrs = {
        "output": attr.string(
            mandatory = True,
        ),
        "generator": attr.label(
            mandatory = True,
            executable = True,
            cfg = "exec",
        ),
    },
)
