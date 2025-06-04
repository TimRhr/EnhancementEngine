# Frequently Asked Questions

## Do I need an email address to use the engine?
Yes. An email address is required for NCBI requests and must be supplied when creating `EnhancementEngine` or running the CLI.

## Can I run analyses from the command line?
Absolutely. The `enhancement-engine` command mirrors the Python API. Use `enhancement-engine --help` for a list of commands.

## Where are results stored?
By default results are kept in memory. Sequence data may be cached under `data/cached_sequences/`. You can clear the cache using `enhancement-engine cache clear`.

## Is this software suitable for clinical use?
No. The Enhancement Engine is intended for research and educational purposes only and should not be used for real medical decisions.
