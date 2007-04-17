####################################################
# override SharedLibrary to only update the libs in '#shlib/' if the symbols
# change.  This is to avoid rebuilding binaries when a shared library has
# changes.
# This would be MUCH simpler to implement if we could override the
# signature to used for the SharedLibrary node itself directly. That is
# certainly possible, but would rely on the internal structure of SCons.

def fasterSharedLibrary(env, library, sources, **args):
	# SCons version compatibility
	if type(library) != type([]):
		library = [library]
        # use the 'quicker' shallow copy method!
        envContentSig=env.Copy()
        envContentSig.TargetSignatures('content')

        cat=env.OriginalSharedLibrary(library, sources)

        # copy all the latest libraries to ONE directory..
        # for our convenience. Could modify the above to
        # build directly to this dir instead.
        catLib = env.Install('#lib', cat) #the CURRENT lib dir

        # now generate the 'interface' file, using the
        # content signature for its target
        catIF=envContentSig.Command(
		'%s.if' % library[0],
                catLib,
                'nm --extern-only $SOURCES | cut -c 12- | sort > $TARGET')

        # install command to copy lib to shlib, where the link
        # actually occurs.  Explicitly make this depend only on
        # the IF file, which has a target content signature.
        # ie only if the Global Symbol list changes, is copied and this the
        # Programs it relinked.
        catLink=env.Command(
		'#linklib/${SHLIBPREFIX}%s${SHLIBSUFFIX}' % library[0],
                '',
                Copy('$TARGET', catLib))

        #Dir('#lib')
        envContentSig.Depends(catLink, catIF)

        #global libs
        #libs += catLib

        return cat

# declaring OriginalSharedLibrary is a bit marginal.  Probably should use
# a functor style object so we can store it in side the object?
#env['BUILDERS']['OriginalSharedLibrary'] = env['BUILDERS']['SharedLibrary']
#env['BUILDERS']['SharedLibrary'] = fasterSharedLibrary
