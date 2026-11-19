"""`python -m idtrack_reproducibility.corpus <subcommand>`: the corpus's steps (pipeline.py)."""
import sys

from idtrack_reproducibility.corpus.pipeline import main

sys.exit(main())
