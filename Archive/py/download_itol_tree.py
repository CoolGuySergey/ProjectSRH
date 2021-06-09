from itolapi import ItolExport


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', help="output file", required=True, type=str)
    parser.add_argument('-t', '--tree', help="tree id", required=True, type=str)
    parser.add_argument('-f', '--format', help="file format", required=True, type=str)
    parser.add_argument('-c', '--column', help="column of interest", required=True, type=str)
    parser.add_argument('--columns', help="available columns", nargs='+', type=str)
    params = parser.parse_args()

    itol_exporter = ItolExport()
    itol_exporter.set_export_param_value('tree', params.tree)
    assert params.format in ['png', 'svg', 'eps', 'ps', 'pdf', 'nexus', 'newick']
    itol_exporter.set_export_param_value('format', params.format)
    itol_exporter.set_export_param_value('display_mode', 2)
    itol_exporter.set_export_param_value('vertical_shift_factor', 0.09)
    itol_exporter.set_export_param_value('label_display', 0)
    col_id = sorted(params.columns).index(params.column)
    itol_exporter.set_export_param_value('datasets_visible', len(params.columns) + col_id)
    itol_exporter.set_export_param_value('tree_x', 250)
    itol_exporter.set_export_param_value('internal_scale', 1)
    itol_exporter.set_export_param_value('internalScale1', 1000)
    itol_exporter.set_export_param_value('internalScale2', 10)
    itol_exporter.set_export_param_value('arc', 358)
    itol_exporter.set_export_param_value('line_width', 2)

    itol_exporter.export(params.output)

