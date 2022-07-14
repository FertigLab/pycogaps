'''
This file runs PyCoGAPS using user-specified parameter inputs from params.yml
This is intended for running development version of PyCoGAPS (from github)
'''

if __name__ == '__main__':
    from pycogaps.parameters import *
    from pycogaps.pycogaps_main import CoGAPS

    import yaml
    import pickle
    print("This vignette was built using pycogaps version", getVersion())

    # read parameter file
    with open("params.yaml", "r") as file:
        prm = yaml.safe_load(file)
    
    # if using AWS bucket server
    aws_prm = prm['aws_params']
    if aws_prm['useAWS']:
        import boto3
        s3 = boto3.client('s3')
        with open(prm['path'], 'wb') as f:
            s3.download_fileobj(aws_prm['downloadBucket'], aws_prm['downloadKey'], f)
    
    # create CoParams object
    params = CoParams(path=prm['path'], transposeData=prm['run_params']['transposeData'], 
                      hdfKey=prm['additional_params']['hdfKey'], hdfRowKey=prm['additional_params']['hdfRowKey'],
                      hdfColKey=prm['additional_params']['hdfColKey'])
    
    # set all standard, sparsity, additional parameters
    setParams(params, prm['standard_params'])
    setParams(params, prm['sparsity_params'])
    setParams(params, prm['additional_params'])

    # set fixed patterns from additional params
    if prm['additional_params']['fixedPatterns'] is not None:
        params.setFixedPatterns(fixedPatterns=prm['additional_params']['fixedPatterns'], whichMatrixFixed=prm['additional_params']['whichMatrixFixed'])

    # set distributed parameters
    dist_prm = prm['distributed_params']
    setParam(params, 'distributed', dist_prm['distributed'])
    if dist_prm['distributed'] is not None:
        params.setAnnotationWeights(annotation=dist_prm['samplingAnnotation'], weight=dist_prm['samplingWeight'])
        params.setDistributedParams(nSets=dist_prm['nSets'], cut=dist_prm['cut'], minNS=dist_prm['minNS'], maxNS=dist_prm['maxNS'])

    # run CoGAPS
    result = CoGAPS(prm['path'], params)

    # save CoGAPS result
    print("Pickling...", end='\r')
    pickle.dump(result, open(prm['result_file'], "wb"))
    print("Pickling complete!")

    if aws_prm['useAWS']:
        with open(prm['result_file'], 'rb') as data:
            s3.upload_fileobj(data, aws_prm['uploadBucket'], prm['uploadKey'])
    