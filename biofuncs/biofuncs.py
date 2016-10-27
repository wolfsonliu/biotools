import pandas as pd
import numpy as np

def gfftable(filename, fileformat, columns):
    '''Parse gff3 or gtf file into pandas DataFrame.'''
    # {{{

    if fileformat == 'gtf': # make different separation according to fileformat.
        sep1, sep2 = ' "', '";'
    elif fileformat == 'gff':
        sep1, sep2 = '=', ';'
    data = pd.read_table(
        filename,
        sep     = '\t',
        header  = None,
        names   = columns,
        index_col = False,
        comment = '#'
    )                           # read gtf/gff data
    attr_split = data[columns[-1]].apply(
        lambda row: [x.split(sep1)
                     for x in row.split(sep2) if len(x) > 2]
    )                      # split the attributes in attribute columns
    attr_dict = attr_split.apply(
        lambda row: dict(
            zip(
                [x[0].strip() for x in row],
                [x[1].strip() for x in row]
            )
        )
    )                           # make splited attributes into dicts
    attr_columns = attr_dict.apply(
        lambda row: list(row.keys())
    ).tolist()                  # get attr columns names
    attr_names = list(
        set([
            attr_columns[i][j] for i in range(len(attr_columns))
            for j in range(len(attr_columns[i]))
        ])
    )                           # get attr columns names
    attr = pd.DataFrame(
        dict(
            zip(
                attr_names,
                [
                    pd.Series(
                        [x[attr_name] if attr_name in x else np.NaN
                         for x in attr_dict]
                    ) for attr_name in attr_names
                ]
            )
        )
    )                           # make attr columns
    data = data.join(
        attr
    )                           # link attr columns with annotation.
    return data

    # }}}

