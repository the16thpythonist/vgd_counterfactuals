import os
import json
from typing import List

import click
from pycomex.cli import ExperimentCLI

from vgd_counterfactuals.utils import PATH
from vgd_counterfactuals.utils import get_version
from vgd_counterfactuals.utils import CsvString

@click.group()
@click.option('--version', is_flag=True, help='Show version and exit.')
@click.pass_context
def cli(ctx: click.Context, version: bool):
    """Command line interface for VGD Counterfactuals."""
    if version:
        click.echo(f'Version: {get_version()}')
        ctx.exit()

    return


@cli.command()
@click.argument('graph_representation', type=str)
@click.option('--k', default=1, show_default=True, type=int, help='Number of hops for the neighborhood.')
@click.option('--domain', default='molecules', show_default=True, type=click.Choice(['molecules', 'color_graphs']), help='Domain type.')
@click.pass_context
def variants(ctx: click.Context, graph_representation: str, k: int, domain: str):
    """Outputs a list of the k-hop neighbors for the given graph representation."""
    if domain == 'molecules':
        from vgd_counterfactuals.generate.molecules import get_neighborhood, DEFAULT_ATOM_VALENCE_MAP
        neighbors = [
            {'value': graph_representation, 'type': 'original', 'org': (), 'mod': ()}
        ]
        for _ in range(k):
            next_neighbors = []
            for data in neighbors:
                # get_neighborhood returns a list of dicts
                next_neighbors.extend(get_neighborhood(data['value']))
            neighbors = next_neighbors
        # Output only the unique neighbor values
        unique_neighbors = {data['value']: data for data in neighbors}.values()
        click.echo(json.dumps(list(unique_neighbors), indent=2))
    elif domain == 'color_graphs':
        from vgd_counterfactuals.generate.colors import get_neighborhood
        neighbors = [
            {'value': graph_representation, 'type': 'original', 'org': (), 'mod': ()}
        ]
        for _ in range(k):
            next_neighbors = []
            for data in neighbors:
                next_neighbors.extend(get_neighborhood(data['value']))
            neighbors = next_neighbors
        unique_neighbors = {data['value']: data for data in neighbors}.values()
        click.echo(json.dumps(list(unique_neighbors), indent=2))
    else:
        click.echo(f'Unknown domain: {domain}')


if __name__ == '__main__':
    cli()
