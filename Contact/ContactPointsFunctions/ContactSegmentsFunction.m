function [ContactSegments, ContactSegmentsElements] = ContactSegmentsFunction(ContactBody)

    nloc_cont = ContactBody.nloc;
    ContactNode_cont = ContactBody.contact.nodalid(:);

    % Each row contains the global node numbers of one segment.
    ContactSegments = [ContactNode_cont(1:end-1), ContactNode_cont(2:end)];

    numberOfSegments = size(ContactSegments,1);
    ContactSegmentsElements = zeros(numberOfSegments,1);

    for ii = 1:numberOfSegments

        node_a = ContactSegments(ii,1);
        node_b = ContactSegments(ii,2);

        % Find the element containing both segment nodes.
        ElementNumber = find(any(nloc_cont == node_a,2) & any(nloc_cont == node_b,2));

        if numel(ElementNumber) ~= 1
            error('Contact segment [%d, %d] must belong to one element.', node_a, node_b);
        end

        ContactSegmentsElements(ii) = ElementNumber;

    end
end