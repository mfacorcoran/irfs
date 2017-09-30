def generate(env, **kw):
    env.Tool('irfInterfaceLib')
    env.Tool('irfUtilLib')
    env.Tool('irfLoaderLib')
    env.Tool('testResponseLib')
    if baseEnv.get('CONTAINERNAME','') != 'ScienceTools_User':
        env.Tool('dc1aResponseLib')
        env.Tool('handoff_responseLib')

def exists(env):
    return 1
