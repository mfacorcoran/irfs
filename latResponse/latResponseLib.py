#$Id$
def generate(env, **kw):
    env.Tool('addLibrary', library=['latResponse'])
    env.Tool('irfInterfaceLib')
    env.Tool('astroLib')
    env.Tool('st_facilitiesLib')
    env.Tool('st_streamLib')
    env.Tool('irfUtilLib')
    env.Tool('tipLib')
    env.Tool('addLibrary', library=env['clhepLibs'])

def exists(env):
    return 1
