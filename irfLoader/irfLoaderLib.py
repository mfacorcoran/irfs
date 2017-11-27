#$Id$
def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['irfLoader'])
    env.Tool('st_streamLib')
    env.Tool('tipLib')
    env.Tool('facilitiesLib')
    env.Tool('st_facilitiesLib')
    env.Tool('irfInterfaceLib')
    env.Tool('latResponseLib')
    env.Tool('testResponseLib')
    if env.get('CONTAINERNAME', '') != 'ScienceTools_User':
        env.Tool('handoff_responseLib')
        env.Tool('dc1aResponseLib')

def exists(env):
    return 1
