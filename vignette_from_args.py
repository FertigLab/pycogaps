'''
This file runs PyCoGAPS using user-specified parameter inputs from params.yml
This is intended for running Docker image of PyCoGAPS (latest pushed)  
'''

if __name__ == '__main__':
    from PyCoGAPS.config import *
    from PyCoGAPS.parameters import *
    from PyCoGAPS.pycogaps_main import CoGAPS

    import yaml
    import pickle
    print("This vignette was built using pycogaps version", getVersion())

    # get parameter file from command line input
    params_file = sys.argv[1]
    PWD = '/'.join(params_file.split('/')[:-1]) + '/'   # leveraging $PWD tag in the run line
    outdir = PWD + 'output/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # read parameter file
    with open(params_file, "r") as file:
        prm = yaml.safe_load(file)
    
    # if using AWS bucket server
    aws_prm = prm['aws_params']
    if aws_prm['useAWS']:
        import boto3
        s3 = boto3.client('s3')
        with open(prm['path'], 'wb') as f:
            s3.download_fileobj(aws_prm['downloadBucket'], aws_prm['downloadKey'], f)
    
    # create CoParams object

    # Note: since data_path=PWD+prm['path'], the path supplied in param.yaml must be relative to the working directory
    data_path = PWD+prm['path']
    
    params = CoParams(path=data_path, transposeData=prm['run_params']['transposeData'], 
                      hdfKey=prm['additional_params']['hdfKey'], hdfRowKey=prm['additional_params']['hdfRowKey'],
                      hdfColKey=prm['additional_params']['hdfColKey'])
    
    # set all standard, sparsity, additional parameters
    setParams(params, prm['standard_params'])
    setParams(params, prm['run_params'])
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
    result = CoGAPS(data_path, params)

    # save CoGAPS result
    print("Pickling...", end='\r')
    # pickle.dump(result, open(prm['result_file'], "wb"))
    pickle.dump(result, open(outdir+prm['result_file'], "wb"))
    print("Pickling complete!")

    if aws_prm['useAWS']:
        with open(prm['result_file'], 'rb') as data:
            s3.upload_fileobj(data, aws_prm['uploadBucket'], prm['uploadKey'])
    