from setuptools import Extension, setup

setup(
    name = "wbgt",
    ext_modules=[
        Extension(
            name="wbgt.ext",  # as it would be imported
                                  # may include packages/namespaces separated by `.`
            sources=["src/wbgt.c"], # all sources are compiled into a single binary file
        ),
    ],
    packages=['wbgt'],
    author = "Mahe Perrette",
    author_email = "mahe.perrette@gmail.com",
    license = "MIT",
)